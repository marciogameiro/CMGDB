import csv

def LoadMorseSetFile(morse_set_fname):
    morse_sets = []
    with open(morse_set_fname, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, quoting=csv.QUOTE_NONNUMERIC, skipinitialspace=True)
        for row in csv_reader:
            morse_sets.append(row)
    return morse_sets
