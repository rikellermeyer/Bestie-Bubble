#!/usr/bin/env python3

import sys

file_input = sys.argv[1]


color_dict = {}
animal_dict = {}
food_dict = {}
data_list = [color_dict, animal_dict, food_dict]
with open(file_input, 'r') as dataset:
    for line in dataset:
        line = line.rstrip()
        [name, color, animal, food]= line.split('\t')
        color_dict[name] = color
        animal_dict[name] = animal
        food_dict[name] = food

name_list = list(color_dict.keys())

user_data = {'color':'red', 'animal':'dog', 'food':'pizza'}

for name in name_list:
    set_overlap = 
    jaccard = set_overlap / set_sum



        
        
    




