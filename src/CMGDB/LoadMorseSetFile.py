import csv

def LoadMorseSetFile(morse_set_fname):
    morse_sets = []
    with open(morse_set_fname, 'r') as csvfile:
        csvreader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC, skipinitialspace=True)
        headers = next(csvreader)
        for row in csvreader:
            morse_sets.append(row)
    return morse_sets
