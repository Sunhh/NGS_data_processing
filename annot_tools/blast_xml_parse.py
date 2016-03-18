#!/usr/bin/python

# import statments
import re
import os
import sys
import argparse
#import time

# Copyright(C) 2013 by David Ream
# Distributed under the Biopython licence. (http://www.biopython.org/DIST/LICENSE) 
# Do not remove this comment.
# if you used this cite the paper, it has some location... which i will add later.

def make_query_dict(query):
    result = {}
    cnt = 1
    for line in open(query).readlines():
        if line[0] == '>':
            name = line[1:].split(' ')[0]
            index = str(cnt)
            result.update({index:name})
            cnt = cnt + 1
    return result

# The query file should be the one that you used to query the database. It HAS to be in fasta format.
# Any other type of format will cause horrible errors.  You should have used a fasta file anyway, so this is
# unlikely to be an issue.
def parse_ncbi_xml_to_csv(infile, outfile, query, delim = "\t"):
    # Bunch of compiled regular expressions.  While there are faster ways to do this, this way is more resistent to  
    # future changes in the BLAST XML output.
    
    # Calling this function will create a dictionary that will link the number of the query with the sequence name
    # so that a more informative query name than say, 'Query1' will be available.
    query_dict = {}
    if query != '':
        query_dict = make_query_dict(query)
    
    result = []
    header_list = ['query id',  'subject id', '% identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score', 'subject description']
    header = delim.join(header_list)
    result.append(header)
    
    # In the statements below, I compile a regular expression. The name of each variable is corresponds to the XML, and is 
    # descriptive of the datum contained. hsp = high scoring pair. 
    
    r_query_id = re.compile("<Iteration_query-ID>")
    r_query_def = re.compile("<Iteration_query-def>")
    r_query_len = re.compile("<Iteration_query-len>")
    
    r_hit_num = re.compile("<Hit_num>")
    r_hit_id = re.compile("<Hit_id>")
    r_hit_def = re.compile("<Hit_def>")
    r_hit_accession = re.compile("<Hit_accession>")
    r_hit_len = re.compile("<Hit_len>")
    
    r_hsp_num = re.compile("<Hsp_num>")
    r_hsp_bit_score = re.compile("<Hsp_bit-score>")
    r_hsp_score = re.compile("<Hsp_score>")
    r_hsp_evalue = re.compile("<Hsp_evalue>")
    r_hsp_query_from = re.compile("<Hsp_query-from>")
    r_hsp_query_to = re.compile("<Hsp_query-to>")
    r_hsp_hit_from = re.compile("<Hsp_hit-from>")
    r_hsp_hit_to = re.compile("<Hsp_hit-to>")
    r_hsp_query_frame = re.compile("<Hsp_query-frame>")
    r_hsp_hit_frame = re.compile("<Hsp_hit-frame>")
    r_hsp_identity = re.compile("<Hsp_identity>")
    r_hsp_positive = re.compile("<Hsp_positive>")
    r_hsp_gaps = re.compile("<Hsp_gaps>")
    r_hsp_align_len = re.compile("<Hsp_align-len>")
    # These variables are unlikely to be useful in a csv due to their huge length.
    r_hsp_qseq = re.compile("<Hsp_qseq>")
    r_hsp_hseq = re.compile("<Hsp_hseq>")
    r_hsp_midline = re.compile("<Hsp_midline>")

    # This is a really boring and long way to split all of the information. I am sorry.
    # doing things this way hopefully makes it easier to understand how to modify my code
    # for your own needs.
    for line in open(infile).readlines():
        if r_query_id.search(line) is not None:
            query_id =  line.split("<Iteration_query-ID>")[1].split("<")[0]
            if query != '':
                query_id = query_dict[query_id.split('_')[1]]
        elif r_query_def.search(line) is not None:
            query_def = line.split("<Iteration_query-def>")[1].split("<")[0]
        elif r_query_len.search(line) is not None:
            query_len = line.split("<Iteration_query-len>")[1].split("<")[0]
        elif r_hit_num.search(line) is not None:
            hit_num = line.split("<Hit_num>")[1].split("<")[0]
        elif r_hit_id.search(line) is not None:
            hit_id = line.split("<Hit_id>")[1].split("<")[0]
        elif r_hit_def.search(line) is not None:
            hit_def = line.split("<Hit_def>")[1].split("<")[0] 
   
        elif r_hit_accession.search(line) is not None:
            hit_accession = line.split("<Hit_accession>")[1].split("<")[0]
        elif r_hit_len.search(line) is not None:
            hit_len = line.split("<Hit_len>")[1].split("<")[0]
        elif r_hsp_num.search(line) is not None:
            hsp_num = line.split("<Hsp_num>")[1].split("<")[0]
        elif r_hsp_bit_score.search(line) is not None:
            hsp_bit_score = line.split("<Hsp_bit-score>")[1].split("<")[0]
            
        elif r_hsp_score.search(line) is not None:
            hsp_score = line.split("<Hsp_score>")[1].split("<")[0]
        elif r_hsp_evalue.search(line) is not None:
            hsp_evalue = line.split("<Hsp_evalue>")[1].split("<")[0]
        elif r_hsp_query_from.search(line) is not None:
            hsp_query_from = line.split("<Hsp_query-from>")[1].split("<")[0]
        elif r_hsp_query_to.search(line) is not None:
            hsp_query_to = line.split("<Hsp_query-to>")[1].split("<")[0]
        elif r_hsp_hit_from.search(line) is not None:
            hsp_hit_from = line.split("<Hsp_hit-from>")[1].split("<")[0]
            
        elif r_hsp_hit_to.search(line) is not None:
            hsp_hit_to = line.split("<Hsp_hit-to>")[1].split("<")[0]
        elif r_hsp_query_frame.search(line) is not None:
            hsp_query_frame = line.split("<Hsp_query-frame>")[1].split("<")[0]
        elif r_hsp_hit_frame.search(line) is not None:
            hsp_hit_frame = line.split("<Hsp_hit-frame>")[1].split("<")[0]
            
        elif r_hsp_identity.search(line) is not None:
            hsp_identity = line.split("<Hsp_identity>")[1].split("<")[0]
        elif r_hsp_positive.search(line) is not None:
            hsp_positive = line.split("<Hsp_positive>")[1].split("<")[0]
        elif r_hsp_gaps.search(line) is not None:
            hsp_gaps = line.split("<Hsp_gaps>")[1].split("<")[0]
        elif r_hsp_align_len.search(line) is not None:
            hsp_align_len = line.split("<Hsp_align-len>")[1].split("<")[0]
            
        elif r_hsp_qseq.search(line) is not None:
            hsp_qseq = line.split("<Hsp_qseq>")[1].split("<")[0]
        elif r_hsp_hseq.search(line) is not None:
            hsp_hseq = line.split("<Hsp_hseq>")[1].split("<")[0]
        elif r_hsp_midline.search(line) is not None:
            hsp_midline = line.split("<Hsp_midline>")[1].split("<")[0]
            # Once we have hit this tag, we have all the variables we will have for the hit.
            # At this point i calculate a few things, pretty up their output, and output this to a file.
            
            # Note: From what I can tell on the web, BLAST calculates PID by dividing the number of identities by the length of the aligned region.
            # If this is wrong for what you are doing modify this line appropriately.
            percent_ident = "%5.2f" %  (((float(hsp_identity)/float(hsp_align_len))*100))
            mismatch = str(int(hsp_align_len) - int(hsp_identity))
            # So hit_def can contain commas, which screw up the csv output. I replace them with colons because this is how I ride.
            hit_def = hit_def.replace(',', ':')
                       
            ##################################### HOW TO CUSTOMIZE YOUR OUTPUT!!!!###########################################################
            # You notice the variable var_list. This is a list, which contains the variables that are interesting                           #
            # for retention. If you do not like the ones I have chosen, please comment the line out, and make your own.                     #
            # To do this: var_list = []                                                                                                     #
            # inside the '[]' brackets, put the name of the variables you want, in the order that you want them to be in.                   #
            # comma seperate the variables in the brackets.                                                                                 #
            # ex. var_list = [var1, var2, var3]  Then you are done!  The header line (where the variables are named) will be wrong though.  #
            # just look above in the code for header_list.  Name each field in the convention provided. Avoiding spaces is advised, but     #
            # obviously I ignored that. Just know that when you open the csv file you should uncheck space as a delimiter.                  #
            #################################################################################################################################
            
            var_list = [query_def, hit_id, percent_ident, hsp_align_len, mismatch, hsp_gaps, hsp_query_from, hsp_query_to, hsp_hit_from, hsp_hit_to, hsp_evalue, hsp_bit_score, hit_def]
            result.append(delim.join(var_list))
    handle = open(outfile, 'w')
    handle.write('\n'.join(result))
    handle.close()
            
# this is quick and dirty, so sorry about how horribly this piece of code was written.

def report_best_hit_on_query(infile, outfile = './best_hit.csv'):
    result = []
    
    # if you are planning on modifying this code:
    # The header list simply is a list of fields that you would like to have in the csv report. 
    # Keep in mind here that the names for the fields cannot have a comma.
    header_list = ['query id',  'subject id', '% identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score', 'subject description']
    header = ','.join(header_list)
    result.append(header)
    
    cur_query = ''
    for line in [i.strip() for i in open(infile).readlines()[1:]]:
        query = line.split(',')[0]
        # Basically this works because the output is ordered by the best hit appearing first.
        # So i simply take the first occurance of any new query hit, and is the best hit.
        if query != cur_query:
            result.append(line)
            cur_query = query
            
    handle = open(outfile, 'w')
    handle.write('\n'.join(result))
    handle.close()

# Check the input values to make sure that they the proper file type.
# For now, it only checks out the file extension, and if the input files exist.
def input_checker(parsed_args):
    infile = parsed_args.infile
    
    # If the input file does not exist provide an error message and exit.
    if not os.path.isfile(parsed_args.infile):
        print "The input file %s does not exist, the program is now exiting." % parsed_args.infile
        sys.exit()
    
    # Check the input file's extension. If not .xml the input is probably the wrong file. 
    # Print an error and exit.    
    f_name, f_extension = os.path.splitext(parsed_args.infile)
    if f_extension != '.xml':
        print "The input file does not have a .xml file extension. The input must be in XML format."
        sys.exit()
        
    # Check that if a query file has been provided, it exists. I am not checking the file extension at this time.
    # People take a lot of liberties with file extensions on fasta files in my experience.    
    query = parsed_args.query
    
    if query != '':
        if not os.path.isfile(parsed_args.infile):
            print "The query file %s was specified, but does not exist, the program is now exiting." % query
            sys.exit()
    
    # Bit of code to auto-generate an output file name based on the checked input.
    # While technically not something to check, I put this code here anyway so that 
    # all option processing is done in one place.
    if parsed_args.outfile == '':
        outfile = infile.split('.')[-2] + '.csv'
    else:
        outfile = parsed_args.outfile
    
    return infile, outfile, query
            
def main():
    # I time the program to tweak the implementation so I know the choices I made run faster.
    # start = time.time()

    parser = argparse.ArgumentParser(description='Convert a BLAST XML file into a CSV for further analysis.')

    parser.add_argument("-i", "--infile", required=True, dest="infile",
                help="Input BLAST XML file. This option is required", metavar="FILE")
                
    parser.add_argument("-o", "--outfile", dest="outfile",
                help="Store the result of the program, this file name should end in '.csv'. If omitted the program will use the base file name and '.csv will be added.",
                metavar="FILE", default='')
                
    parser.add_argument("-q", "--query", dest="query", metavar="FILE", default='',
                help="Query fasta file previously used by BLAST to create the input XML file. This file is optional, but will make the program output more readable.")

    infile, outfile, query = input_checker(parser.parse_args())
    
    # This is where the computational component of the program begins.
    # parse_ncbi_xml_to_csv() and report_best_hit_on_query() do all the work.
    parse_ncbi_xml_to_csv(infile, outfile, query)

    best_hit_infile = outfile
    best_hit_outfile = outfile.replace('.', '_best_hits.')
    
    report_best_hit_on_query(best_hit_infile, best_hit_outfile)
    
    #end = time.time()
    #print end - start

if __name__ == '__main__':
    main()
