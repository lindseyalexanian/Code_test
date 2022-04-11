#!/usr/bin/env python3
# biointech.py


"""This script counts the number of high-scoring peptides in each gene and output the scaled (per unit length) 
and ranked values for each gene in a csv file."""


import pandas as pd
import numpy as np
import csv

def get_header_and_sequence_lists(file_handle):
    """This function separates headers and sequences into separate lists, and also
    ensures that there is an equal number of each"""
    header_lines = []
    sequence_lines = []
    seq = ''
    with open(file_handle, "r") as f:

        try:
            lines = f.readlines()
            for line in lines:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if len(header_lines) > 0 and seq != '':
                        sequence_lines.append(seq)
                    seq = ''
                    header_lines.append(line)

                else:
                    seq += line

            seq2 = seq.rstrip('\n')
            if seq2:
                sequence_lines.append(seq2)

        except AttributeError:
            for i in file_handle:
                i = i.rstrip('\n')
                if i.startswith('>'):
                    if len(header_lines) > 0 and seq != '':
                        sequence_lines.append(seq)
                    seq = ''
                    header_lines.append(i)

                else:
                    seq += i

            seq2 = seq.rstrip('\n')
            if seq2:
                sequence_lines.append(seq2)

        
        return header_lines, sequence_lines

def peps(filename):
    """This function extracts the peptide dataframe, first column,
    and second column. It also makes a dictionary w/ the first column as key and 
    second column as value. It returns all four of these."""

    # read file
    peptides = pd.read_csv(filename, header=0)

    # first column
    first = peptides['peptide'].tolist()

    # second column
    second = peptides['score'].tolist()

    # make a dictionary
    csv_to_dict = {first[i]: second[i] for i in range(len(first))}


    return peptides, first, second, csv_to_dict

class Pep_stats:
    def __init__(self, column):
        self.column = column

    def median_pep(self):
        medianpep = np.median(self.column)
        return medianpep

    def mean_pep(self):
        meanpep = np.mean(self.column)
        return meanpep

    def min_pep(self):
        minimum = min(self.column)
        return minimum

    def max_pep(self):
        maximum = max(self.column)
        return maximum

    def percentile(self):
        per_75 = np.percentile(self.column, 75)
        return per_75


def get_counts_of_pep():

    peptides, first, second, pep_dict = peps("peptide-scores.csv")

    headers, seq = get_header_and_sequence_lists("genes.fa")

    # initialze lst for count of peptides in each sequence
    num_of_pep = []
    for item in seq:
        # initialize couter
        counter = 0
        for pep in first:
            if pep in item:
                if item.count(pep) > 1:
                    counter += item.count(pep)
                else:
                    counter += 1
        
        # + count to list
        num_of_pep.append(counter)

    # combine the sequences names (headers) and the peptide counts into a dictionary
    pep_count_dict = {headers[i]: num_of_pep[i] for i in range(len(headers))}

    # the function returns the dictionary
    return pep_count_dict

def _high_scoring(dictionary, column):
    # initialize list for only high scores
    high_scoring_peps =[]

    # get the 75th percentile
    statistics = Pep_stats(column)
    per75 = statistics.percentile()

    for item in column:
        if item >= per75:
            high_scoring_peps.append(item)
    
    return high_scoring_peps


def high_scoring(dictionary, column):
    # get the list of scores
    high_scoring_pep_scores = _high_scoring(dictionary, column)

    # initialize list for the high scoring peps
    pep_names_75 = []
        
    # loop through the dictionary of peps to see which ones have the high scoring ones
    for key in dictionary:
        if dictionary[key] in high_scoring_pep_scores:
            pep_names_75.append(key)

    # return list of pep names
    return pep_names_75



def get_counts_high_scoring_pep(high_scoring_peps_only_list, sequence, headers):
    # initialize list for count of pep in each
    num_of_pep2 = []
    for item in sequence:
        # initialize counter
        counter = 0
        for pep in high_scoring_peps_only_list:
            if pep in item:
                if item.count(pep) > 1:
                    counter += item.count(pep)
                else:
                    counter += 1
        
        # add count to list
        num_of_pep2.append(counter)

    # combine the sequences names (headers) and the pep counts into a dictionary
    pep_count_dict2 = {headers[i]: num_of_pep2[i] for i in range(len(headers))}

    return pep_count_dict2


def lets_scale_it(pep_count_dict2, sequence, headers):
    """Scales"""
    # initialize seq length list
    seq_length = []
    for item in sequence:
        length_of_seq = len(item)
        seq_length.append(length_of_seq)

    #print(seq_length)

    # new scaled dict
    extract_value = list(pep_count_dict2.values())

    #print(extract_value)

    scaled_vals = [i / j for i, j in zip(extract_value, seq_length)]


    # combine the sequences names (headers) and the pep counts into a dictionary
    pep_count_dict_scaled = {headers[i]: scaled_vals[i] for i in range(len(headers))}
    return pep_count_dict_scaled
    

def sort_your_dictionary(dictionary):
    """sort dictionary by value"""
    sorted_list = sorted(dictionary.items(), key=lambda x: x[1], reverse=True)
    sorted_dict = {}
    for item in sorted_list:
        sorted_dict[item[0]] = item[1]

    key_pop_list = []
    for key, val in sorted_dict.items():
        if val == 0:
            key_pop_list.append(key)

    for key in key_pop_list:
        sorted_dict.pop(key)

    return sorted_dict

def make_them_ranks(ordered_list):
    ranks = []
    counter = 0
    for i in range(0, len(ordered_list)):
        
        if i == 0:
            counter += 1
            ranks.append(counter)


        elif ordered_list[i] < ordered_list[i - 1]:
            counter += 1
            ranks.append(counter)
        
        else:
            ranks.append(counter)
        
    
    return ranks



def compile():
    with open("outfile.txt", "w") as outfile:

        headers, seq = get_header_and_sequence_lists("genes.fa")
        outfile.write(f"The number of genes is {len(headers)}" + "\n")

        # call peps function
        df, col1, col2, df_as_dict = peps("peptide-scores.csv")
        outfile.write(f"The peptide count is {len(col1)}" + "\n")

        # write stats to file
        statistics = Pep_stats(col2)
        # store the percentile stat, will be used later
        per75 = statistics.percentile()
        # write it to file
        outfile.write(f"The 75th percentile is at score {per75}" + "\n")
        print(f"The 75th percentile is at score {per75}" + "\n")


        stats_dict = {"Minimum:": statistics.min_pep(), "Maximum:": statistics.max_pep(), 
                    "Median:": statistics.median_pep(), "Mean:": statistics.mean_pep()}
        # print and write each stat to file
        for key in stats_dict:
            print(key, stats_dict[key])
            outfile.write(f"{key}\t{stats_dict[key]}" + "\n")

        # now get pep count for each entry
        pep_dictionary = get_counts_of_pep()
        # print
        print(f"\nBelow are the pep counts.")
        outfile.write(f"\nBelow are the pep counts.\n")
        # write info to outfile
        for key in pep_dictionary:
            print(f"{key[1:]}: {pep_dictionary[key]}")
            outfile.write(f"{key[1:]}:\t{pep_dictionary[key]}" + "\n")

        # get list of pep names
        pep_names_75 = high_scoring(df_as_dict, col2)

        # get higher counts
        pep_count_dict2 = get_counts_high_scoring_pep(pep_names_75, seq, headers)
        print(f"\nBelow are the 75th per. + counts.")
        outfile.write(f"\nBelow are the 75th per. + counts.\n")
        for key in pep_count_dict2:
            print(f"{key[1:]}: {pep_count_dict2[key]}")
            outfile.write(f"{key[1:]}:\t{pep_count_dict2[key]}" + "\n")


        scaled_dict = lets_scale_it(pep_count_dict2, seq, headers)
        #print(scaled_dict)
        sorted_scaled_dict = sort_your_dictionary(scaled_dict)
        #print(sorted_scaled_dict)
        #print(rank_nums)

        print(f"\nBelow are the scales!")
        outfile.write(f"\nBelow are the scales!\n")
        for key in sorted_scaled_dict:
            print(f"{key[1:]}: {sorted_scaled_dict[key]}")
            outfile.write(f"{key[1:]}:\t{sorted_scaled_dict[key]}" + "\n")


        ordered_headers = list(sorted_scaled_dict.keys())
        ordered_scales = list(sorted_scaled_dict.values())

        # get the ranks in a list
        rank_nums = make_them_ranks(ordered_scales)


        # make dict with ranks (name and rank)
        dict_scaled_ranked = {ordered_headers[i]: rank_nums[i] for i in range(len(ordered_headers))}

        # print and write to file
        print(f"\nBelow are the ranks!")
        outfile.write(f"\nBelow are the ranks!\n")
        for key in dict_scaled_ranked:
            print(f"{key[1:]}: {dict_scaled_ranked[key]}")
            outfile.write(f"{key[1:]}:\t{dict_scaled_ranked[key]}" + "\n")

        return dict_scaled_ranked, sorted_scaled_dict


def main():
    dict_scaled_ranked, sorted_scaled_dict = compile()

    ranked_list_of_lists = []
    for key, val in dict_scaled_ranked.items():
        ranked_list_of_lists.append([key[1:]] + [val])

    print(ranked_list_of_lists)

    csv_header = ["name", "rank"]
    with open('rank.csv', 'w', encoding='UTF8', newline='') as csv_f:
        writer = csv.writer(csv_f)
        
        # write the header
        writer.writerow(csv_header)

        for i in ranked_list_of_lists:
            writer.writerow(i)



    sort_scale_list_of_lists = []
    for key, val in sorted_scaled_dict.items():
        sort_scale_list_of_lists.append([key[1:]] + [val])

    print(sort_scale_list_of_lists)

    csv_header2 = ["name", "scale_unit"]
    with open('scaled.csv', 'w', encoding='UTF8', newline='') as csv_f2:
        writer2 = csv.writer(csv_f2)
        
        # write the header
        writer2.writerow(csv_header2)

        for i in sort_scale_list_of_lists:
            writer2.writerow(i)





# call main function
if __name__ == "__main__":
    main()