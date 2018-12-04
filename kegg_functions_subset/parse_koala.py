#!//anaconda/bin/python

'''

Parses the output file from GhostKoala to retrieve the higher order function (B)
and the number of gene products that are associated to the specific function (D). Removes
duplicate entries present in B 

	Input: ghostkoala BRITE output (.kegg) in .txt format
	Output: list of gene products for each higher-order function (.txt?)

'''

import os
import pandas as pd 

def parse_koala_file(input_f, input_dir):
	"""
	Fill
	"""
	strain_name = input_f[:-10]								#keeping the GCA_XXXXXXX and removing _koala.txt

	#making intermediate output file
	outfile = open("%soutfile_koala_parse.txt" % input_dir, "w")

	#removing the text at the top of of the textfile
	with open ("%s%s" % (input_dir, input_f), "r") as infile:
		for pos, line in enumerate(infile):
			if pos >2:
				if line.startswith("B  "):
					outfile.write(line)
				elif line.startswith("D  "):
					outfile.write(line)
			else:
				continue

	outfile.close()

	#re-reading the intermediate output file and now keeping anything starting with a B or a D
	#storing the whole file to a list for easy access

	all_names = []
	with open("%soutfile_koala_parse.txt" % input_dir, "r") as infile2:
		for line in infile2:
			line = line.replace("\n", "")
			all_names.append(line)

	#counter to track the last b position for the special case next to the break in the while loop below
	b_store = []
	b_count = 0
	for line in all_names:
		if line.startswith("B  "):
			b_count = b_count + 1
			b_store.append(line)

	#storage objects for the while loop
	b_titles = []
	d_titles = []
	temp_titles = []

	i = 0 
	while i < len(all_names):

		if all_names[i].startswith("B  "):
			entry = all_names[i]
			b_titles.append(entry)

			if i != 0:

				d_titles.append(temp_titles)

			temp_titles = []

		else:
			if all_names[i].startswith("D  "):
				entry = all_names[i]
				temp_titles.append(entry)

		if len(b_titles) == b_count:			#the last b title needs to have everything captured after its i position

			temp_titles = all_names[i+1 : ]
			d_titles.append(temp_titles)

			break

		i += 1


	#removing duplicates from each classification level (but kept in different classifications)
	d_titles_cleaned = []
	for i in d_titles:
		new_list = set(i)
		d_titles_cleaned.append(new_list)

	#making a count list of the cleaned gene list
	gene_counts = []
	for genes in d_titles_cleaned:
		gene_counts.append(len(genes))
	
	#turning the classification terms and the counts of genes within into a dataframe
	df = pd.DataFrame(
		{"Categories" : b_titles,
		"Count" : gene_counts,
		"Sample" : strain_name,
		})

	print(strain_name)
	print ("categories in genome:  ", len(df))
	return df

"""

Script begins here 

"""
#input directory
koala_dir = "/<path to 'kegg_functions_subset'>/"
#output directory
out_koala_dir = "/<path to 'kegg_functions_subset'>/parse_output"
#output filename
out_filename = "/<path to 'kegg_functions_subset'>/output/parse_c_level"

#final storage vector
df_store = pd.DataFrame()

#parsing individual koala files and adding the df to the final df_store
for file in os.listdir(koala_dir):
	if file.endswith("koala.txt"):
		temp_df = parse_koala_file(file, koala_dir)

		df_store = pd.concat([df_store, temp_df], ignore_index = True)

	
#writing master df_store to csv
df_store.to_csv("%s%s.csv" % (out_koala_dir, out_filename), index = False)

