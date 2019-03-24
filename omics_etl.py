import threading
import json
from pyensembl import EnsemblRelease
import cx_Oracle
from queue import Queue
import csv
from threading import Thread
import boto3
import datetime
import requests
from pysam import VariantFile
import pysam
from datetime import datetime
from elasticsearch import Elasticsearch, helpers
import logging
import mygene
import sys

es_log = logging.getLogger("elasticsearch")
es_log.setLevel(logging.CRITICAL)
boto_log = logging.getLogger("botocore")
boto_log.setLevel(logging.CRITICAL)


class VCFFileObject:
    def __init__(self, file_location):

        self.file_location = file_location
        self.index_file = (self.file_location + '.tbi')
        self.GenomeBuild = None
        self.info_dict = {}
        self.format_dict = {}
        self.vcf = None
        self.header = None
        self.header_exists = True

    def process_files(self):
        try:
            self.vcf = VariantFile(self.file_location)  # auto-detect input format
            self.header = self.vcf.header
            self.get_genome_build()
            if self.header_exists:
                self.make_info_dict()
                self.make_format_dict()
        except Exception as e:
            print(e, self.file_location)
            pass

    def get_genome_build(self):
        builds = ['hg18', 'hg19', 'hg38', 'human_g1k_v37', 'GRCh38', 'NCBI36']
        for h in self.vcf.header.records:
            h = str(h)
            if h.startswith('##reference'):
                for b in builds:
                    if b in h:
                        self.GenomeBuild = b
                        break
                if self.GenomeBuild == 'human_g1k_v37':
                    self.GenomeBuild = 'hg19'
                if self.GenomeBuild == 'GRCh38':
                    self.GenomeBuild = 'hg38'
                if self.GenomeBuild == "NCBI36":
                    self.GenomeBuild = 'hg18'

        if self.GenomeBuild is None:
            print('No header in file', self.file_location)
            self.header_exists = False
            self.GenomeBuild = None

    def make_info_dict(self):
        for h in self.vcf.header.records:
            h = str(h)
            if h.startswith('##INFO'):
                info_keys = {}
                info = h.split('<')
                info = info[1].split(',')
                ID = info[0].split('=')
                Number = info[1].split('=')
                Type = info[2].split('=')
                Description = info[3].split('=')

                # info_keys[ID[0]] = ID[1]
                info_keys[Number[0]] = Number[1]
                info_keys[Type[0]] = Type[1]
                d = Description[1].split('"')
                info_keys[Description[0]] = d[1]

                self.info_dict[ID[1]] = info_keys

    def make_format_dict(self):
        for h in self.vcf.header.records:
            h = str(h)
            if h.startswith('##FORMAT'):
                format_keys = {}
                format = h.split('<')
                format = format[1].split(',')
                if len(format) > 4:
                    format[3] = ','.join([format[-2], format[-1]])
                    format.remove(format[-1])
                ID = format[0].split('=')
                Number = format[1].split('=')
                Type = format[2].split('=')
                Description = format[3].split('=')

                format_keys[Number[0]] = Number[1]
                format_keys[Type[0]] = Type[1]
                d = Description[1].split('"')
                format_keys[Description[0]] = d[1]

                self.format_dict[ID[1]] = format_keys

class Worker(Thread):
    """ Thread executing tasks from a given tasks queue """
    def __init__(self, tasks):
        Thread.__init__(self)
        self.tasks = tasks
        self.daemon = True
        self.start()

    def run(self):
        while True:
            func, args, kargs = self.tasks.get()
            try:
                func(*args, **kargs)
            except Exception as e:
                # An exception happened in this thread
                print(e)
            finally:
                # Mark this task as done, whether an exception happened or not
                self.tasks.task_done()


class ThreadPool:
    """ Pool of threads consuming tasks from a queue """
    def __init__(self, num_threads):
        self.tasks = Queue(num_threads)
        for _ in range(num_threads):
            Worker(self.tasks)

    def add_task(self, func, *args, **kargs):
        """ Add a task to the queue """
        self.tasks.put((func, args, kargs))

    def map(self, func, args_list):
        """ Add a list of tasks to the queue """
        for args in args_list:
            self.add_task(func, args)

    def wait_completion(self):
        """ Wait for completion of all the tasks in the queue """
        self.tasks.join()

class ProcessOmicsData:
    def __init__(self, VCFFileObject, cytoband, mg, es):
        self.VCFFile = VCFFileObject
        self.cytobandDict = cytoband
        self.mg = mg
        self.es = es
        self.GeneToProteinDict = {}
        self.GenomeBuild = self.VCFFile.GenomeBuild
        self.datahg19 = EnsemblRelease(75)
        self.datahg38 = EnsemblRelease(87)
        self.datahg18 = EnsemblRelease(54)
        self.header = self.VCFFile.header
        self.vcf = self.VCFFile.vcf
        self.info_dict = self.VCFFile.info_dict
        self.format_dict = self.VCFFile.format_dict
        self.organism = "Homo Sapiens"
        self.snpLink = "None"
        self.referenceBuild = 75
        self.server = 'http://rest.ensembl.org'
        self.geneName = "None"
        self.ENSG = "None"
        self.geneLink = "None"
        self.cytoband = "None"
        self.mapping_url = 'http://www.uniprot.org/mapping/'
        self.uniprot_url = 'http://www.uniprot.org/uniprot/'
        self.proteinName = "None"
        self.UniProtLink = "None"
        self.genotype_format = {}
        self.vcf_info = {}
        self.bulk_action = []

    def start_threads(self):
        # Instantiate a thread pool with i worker threads

        pool = ThreadPool(10)

        # Add the jobs in bulk to the thread pool
        startTime = datetime.now()
        geneList = set()
        self.GeneToProteinDict = {}
        for rec in self.vcf.fetch():
            pool.add_task(self.process_omics, rec)
            pool.wait_completion()
            self.geneName = "None"
            self.ENSG = "None"
            self.geneLink = "None"
            self.proteinName = "None"

            # gene and protein extraction functions are not thread safe
            # must be done outside of mulithread function
            self.get_gene_info()
            if self.geneName is not "None":
                geneList.add(self.geneName)
            self.write_omics_file()

            bulk = {
                "_index": 'vcf-ssc',
                "_type": 'vcf-rec',
                "_source": self.omics_data,
            }
            self.bulk_action.append(bulk)

        # this allows for querying proteins in bulk
        self.get_protein_from_gene(geneList)

        # match gene to corresponding protein in bulk_action (list of records)
        self.get_protein_info()

        try:
            helpers.bulk(self.es, self.bulk_action)
            print('Bulk upload successful! Total time to process {} records in {}: {}'.format(len(
                self.bulk_action),
                self.VCFFile.file_location,
                datetime.now() - startTime))

        except Exception as ex:
            print('Error:', ex)


    def process_omics(self, rec):
        self.rec = rec

        self.genotype_format['KEYS'] = self.format_dict
        if self.rec.id is not '.':
            self.snpLink = "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs={}".format(self.rec.id)

        self.chromosome = self.rec.contig
        self.position = self.rec.pos

        self.vcf_info['KEYS'] = self.info_dict
        for key in list(self.rec.info.keys()):
            self.vcf_info[key] = self.rec.info[key]

        chr = ''.join(['chr', self.chromosome])
        for position in self.cytobandDict[self.GenomeBuild]:
            for key, value in position.items():
                if key == chr:
                    for k, v in value.items():
                        if self.position >= int(k):
                            for stop, pos in v.items():
                                if self.position < int(stop):
                                    self.cytoband = pos

    def get_gene_info(self):
        data = self.datahg19
        #if self.GenomeBuild == 'hg38':
         #   data = self.datahg38
        #if self.GenomeBuild == 'hg18':
         #   data = self.datahg18
        try:
            g = data.gene_names_at_locus(contig=self.chromosome, position=self.position)
            if len(g) > 1:
                if '.' in g[0]:
                    self.geneName = g[1]
            else:
                self.geneName = g[0]
            self.ENSG = data.gene_ids_of_gene_name(self.geneName)[0]
            self.geneLink = "http://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g={};r={}".format(self.ENSG, self.chromosome)
        except Exception as e:
            pass


    def get_protein_from_gene(self, geneList):
        global mg

        results = mg.querymany(geneList, scopes='symbol', species=9606, verbose=False)

        for r in results:
            g = r['query']
            try:
               p = r['name']
               self.GeneToProteinDict[g] = p
            except KeyError:
                pass

        #print('Results:', len(results))
        #print('GENE LIST', len(geneList))
        #print('Dict', len(self.GeneToProteinDict))

    def get_protein_info(self):
        for b in self.bulk_action:
            gene = b['_source']['GENE(ENSEMBL)']['geneName']
            if gene is not 'None':
                try:
                    self.protein = self.GeneToProteinDict[gene]
                    b['_source']['PROTEIN(UNIPROT)']['proteinName'] = self.protein
                    #print(gene, self.protein)
                except Exception as e:
                    b['_source']['PROTEIN(UNIPROT)']['proteinName'] = 'None'

    def write_omics_file(self):

        alternate = list(self.rec.alts)
        alternate  = ' '.join(alternate)

        self.omics_data = {"chromosome": self.chromosome,
                           "position": self.position,
                           "referenceAllele": self.rec.ref,
                           "alternateAllele": alternate,
                           "start": self.rec.pos,
                           "end": self.rec.pos,
                           "quality": self.rec.qual,
                           "filter": list(self.rec.filter),
                           "info": self.vcf_info,
                           "genotypeSample": self.genotype_format,
                           "cytobandPosition": self.cytoband,
                           "organism": self.organism,
                           "RSID(DB_SNP)": {
                               "rsid": self.rec.id,
                               "link": self.snpLink,
                           },
                           "GENE(ENSEMBL)": {
                               "geneName": self.geneName,
                               "EnsemblGeneID": self.ENSG,
                               "link": self.geneLink
                           },
                           "PROTEIN(UNIPROT)": {
                               "proteinName": self.proteinName,
                               "link": self.UniProtLink
                           }
                           }


if __name__ == "__main__":

    def create_cytoband():

        print('Pulling cytoband information now...')

        cytobandDict = {}

        # http://hgdownload.cse.ucsc.edu/downloads.html#human
        paths = {'hg19': 'cytoBandhg19.txt',
                 'hg18': 'cytoBandhg18.txt',
                 'hg38': 'cytoBandhg38.txt'}
        for k, v in paths.items():
            path = v
            chrom_list = []

            with open(path, 'r') as tsv_file:
                try:
                    tsv = csv.reader(tsv_file, delimiter="\t")
                    for row in tsv:
                        if '_' not in row[0]:
                            stop = {}
                            start = {}
                            chromosome = {}
                            stop[row[2]] = row[3]
                            start[row[1]] = stop
                            chromosome[row[0]] = start
                            chrom_list.append(chromosome)

                    cytobandDict[k] = chrom_list
                except OSError as e:
                    print('cytoband error', e)

        return cytobandDict

    startTime = datetime.now()

    # Get the service resource
    sqs = boto3.resource('sqs')

    # Create the queue. This returns an SQS.Queue instance
    queue = sqs.create_queue(QueueName='vcf_file_queue', Attributes={'DelaySeconds': '5'})

    file_list = []
    for file in file_list:
        body = json.dumps(file)
        response = queue.send_message(MessageBody=body)
        

    cytoband_info = create_cytoband()
    
    HOST = "" #elasticscearch host 

    es = Elasticsearch([{HOST}], timeout=30, max_retries=10, retry_on_timeout=True)

    mg = mygene.MyGeneInfo()

    messages = queue.receive_messages(MaxNumberOfMessages=10, WaitTimeSeconds=10, VisibilityTimeout=12000)
    while len(messages) > 0:
        for message in messages:
            m = json.loads(message.body)
            message.delete()
            vcf_data = VCFFileObject(m)
            vcf_data.process_files()
            if vcf_data.header_exists:
                omics_data = ProcessOmicsData(vcf_data, cytoband_info, mg, es)
                omics_data.start_threads()
                gp = omics_data.GeneToProteinDict
        messages = queue.receive_messages(VisibilityTimeout=12000)#MaxNumberOfMessages=1, WaitTimeSeconds=10)

    print('Total time: ', datetime.now() - startTime)
