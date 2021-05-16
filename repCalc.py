#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse,subprocess,re

def nucmer(infile,outfile):
    cmd = 'nucmer --nosimplify -c 50 -coords -p {} {} {}'.format(outfile,infile,infile)  #Runs nucmer
    subprocess.call(cmd,shell=True)
    coords(outfile)

def coords(outfile): ###Assign each columm of the .coords file to a defined temporal variable
    fhand = open(outfile + ".coords","r")
    fout = open(outfile + ".repeat.summary","w")
    repeats80 = 0
    repeats90 = 0
    contigs = 0
    for line in fhand:
        if re.search("\s+\d+\s+\d+\s+\|",line):
            line = re.sub("\s+","\t",line)
            line = line.strip().split("|")
            SE1 = line[0].strip().split("\t")
            S1 = SE1[0]
            E1 = SE1[1]
            SE2 = line[1].strip().split("\t")
            S2 = SE2[0]
            E2 = SE2[1]
            LEN = line[2].strip().split("\t")
            LEN1 = LEN[0]
            LEN2 = LEN[1]
            IDEN = line[3]
            TAG = line[4].strip().split("\t")
            TAG1 = TAG[0]
            TAG2 = TAG[1]
            contigs,repeats90,repeats80 = sum(S1,E1,S2,E2,TAG1,TAG2,LEN1,LEN2,IDEN,repeats80,repeats90,contigs)
    frac90 = round((float(repeats90)/float(contigs))*100,2)
    frac80 = round((float(repeats80)/float(contigs))*100,2)
    print  '\nThe percentage of the genome composed by repeats at 90% identity threshold is {}%: {} bp from {} bp.\n'.format(frac90,repeats90,contigs)
    print  '\nThe percentage of the genome composed by repeats at 80% identity threshold is {}%: {} bp from {} bp.\n'.format(frac80,repeats80,contigs)
    fout.write('Genome size\t90%\tfrac90%\t80%\tfrac80%\n')
    fout.write('{}\t{}\t{}\t{}\t{}\n'.format(contigs,repeats90,frac90,repeats80,frac80)) 
    fout.close()  
          
def sum(S1,E1,S2,E2,TAG1,TAG2,LEN1,LEN2,IDEN,repeats80,repeats90,contigs):
    if TAG1 == TAG2 and (S1 == S2 and E1 == E2):
	contigs += int(LEN1)
    else:
         if float(IDEN) >= 90:
             repeats90 += int(LEN1) + int(LEN2)
             repeats80 += int(LEN1) + int(LEN2)
         elif float(IDEN) >= 80:
             repeats80 += int(LEN1) + int(LEN2)   
    return contigs,repeats90,repeats80                          
   
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Use nucmer to calculate the total lenght (bp) spanned by repeated regions.")
    parser.add_argument("-i","--input",action='store',required = True, help = "Input fasta file.")
    parser.add_argument("-o","--output",action='store',required = True, help = "Output basename.")
    args = parser.parse_args()
    infile = args.input
    outfile = args.output
    nucmer(infile,outfile)
