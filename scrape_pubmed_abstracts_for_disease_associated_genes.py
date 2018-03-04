from Bio import Entrez
from Bio import Medline
from textblob import TextBlob
import wikipedia
from wordcloud import WordCloud
import os
import numpy as np
from PIL import Image
import nltk
import sys
from collections import Counter, defaultdict

words = set(nltk.corpus.words.words())
currdir = os.path.dirname(__file__)

def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)

## Initialize Default Parameters ###
MAX_COUNT = 100
DISEASE_TERM = 'DIABETES'

## Read in command line args
if(len(sys.argv) < 2):
   print ('\nRunning with default parameters' )
   MAX_COUNT = 100
   DISEASE_TERM = 'DIABETES'
elif hasNumbers(sys.argv[1]):
   print ('\nUnrecognized Disease Term, Running with default parameters')
   MAX_COUNT = 100
   DISEASE_TERM = 'DIABETES'
elif not hasNumbers(sys.argv[2]):
   print ('\nUnrecognized Disease Term, Running with default parameters')
   MAX_COUNT = 100
   DISEASE_TERM = 'DIABETES'
elif (((len(sys.argv) > 2) and not hasNumbers(sys.argv[1])) and ((int(sys.argv[2])) and (int(sys.argv[2]) < 10000) and (int(sys.argv[2]) > 99))):
   DISEASE_TERM = str(sys.argv[1]).upper()
   MAX_COUNT = int(sys.argv[2])
else:
   print ('\nRunning with default parameters' )
   MAX_COUNT = 100
   DISEASE_TERM = 'DIABETES'

print ('\nMining ', MAX_COUNT, ' Pubmed abstracts for ',DISEASE_TERM,' genes\n' )
TERM = DISEASE_TERM + ' AND GENE'
png_file = DISEASE_TERM + ".png"

def create_wordcloud(text):
        mask = np.array(Image.open(os.path.join(currdir, "template.png")))
        #stopwords = set(STOPWORDS)
        wc = WordCloud(background_color="white",
                                        max_words=1000,
                                        mask=mask,
                                        collocations = False,
                      )
        wc.generate(text)
        wc.to_file(os.path.join(currdir, png_file))
        img = Image.open(png_file)
        img.show()

## Make a dictionary of gene names 
d = {}
with open("FINAL_LIST_PROTEIN_CODING_GENE_HGNC_IDS_unique.txt") as f:
    for line in f:
       new = line.strip()
       for w in nltk.wordpunct_tokenize(new):
          wl = w.lower()
          if wl not in words:
             d[new] = new

# Connect to the biopython database
Entrez.email = 'abhi_r6@hotmal.com'
h = Entrez.esearch(db='pubmed', retmax = MAX_COUNT, sort='relevance', term=TERM)
result = Entrez.read(h)
ids = result['IdList']
h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
records = Medline.parse(h)

str1 = ""
total_parsed_abstracts = 0
total_genes = 0
for record in records:
            pub_abstract = record.get("AB")
            pub_id = record.get("PMID")
            if(pub_abstract):
               analysis=TextBlob(pub_abstract)
               total_parsed_abstracts += 1
               done = {}
               for word in analysis.words:
                 u_word = word.upper()	
                 if u_word in d.keys() and u_word not in done.keys():
                       str1 = str1+" "+u_word
                       total_genes += 1
                 done[u_word] = u_word

print(total_parsed_abstracts, 'Pubmed abstracts are avaialble for mining')

if (total_parsed_abstracts < 10):
   print ('Less than 10 abstracts available, exiting\n')
   exit()

if (total_genes < 5):
   print ('Less than 5 genes, exiting\n')
   exit()

print('\nDescription: ', wikipedia.summary(DISEASE_TERM, sentences=1))

print('\nTop 10 genes found in : ', total_parsed_abstracts, 'Pubmed abstracts\n')

result = defaultdict(lambda: [0, []])
for i, l in enumerate(str1.splitlines()):
    for k, v in Counter(l.split()).items():
        result[k][0] += v
        result[k][1].append(i+1)
items = sorted(result.items(), key=lambda t: (-t[1][0], t[0].lower()))
counter = 1
for (k, v) in items:
    if counter < 11: 
      print('{1} {0}'.format(k, *v))
    counter += 1 

gene_cloud = create_wordcloud(str1)
