#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:03:21 2023

@author: lucas
"""

import os

# Add the path and file name
directory = '/Users/lucas/Documents/UZH/4_jahr_HS_2023/ResearchProject_NMR/Data/Phos_230822_b2_01_GM15_HMQC'
file_name = 'Sample_List.py'
file_path = os.path.join(directory, file_name)  # Create the full path 



# %%

# Function which searches for the keyword and gives the accoarding number + 2 in a list

def process_file_for_keyword(keyword, file_path):
    numbers_list = []

    with open(file_path, 'r') as file:
        for line in file:
            if keyword in line:
                parts = line.strip().split(',')
                if len(parts) >= 3:
                    try:
                        number = int(parts[2].strip())
                        numbers_list.append(number + 2)
                    except ValueError:
                        pass

    return numbers_list


# %%


# Define a list of keywords
keywords = ['b2AR 01 wt', 'b2Ar 02 wt', 'b2AR 03 wt', 'b2AR 01 mut', 'b2AR02 mut', 'b2AR 03 mut']  

# Process the data for each keyword
for keyword in keywords:
    numbers_list = process_file_for_keyword(keyword, file_path)
    numbers_list.reverse()
    print(f"{keyword}, HMQC data", numbers_list)