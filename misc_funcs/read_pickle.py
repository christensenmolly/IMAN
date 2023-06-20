#!/usr/bin/python
import pickle
import argparse


def read_pickle(file):
    with open(file,'rb') as file_object:
        raw_data = file_object.read()

    results = pickle.loads(raw_data)
    return results  



def main(input_file='results.pkl'):
        results = read_pickle(input_file)
        for key, value in results.items(): 
            print('%s: %s' % (key, str(value))) 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parsing a pickle file")
    parser.add_argument("--i", help="Input pickle file", type=str, default='results.pkl')
    args = parser.parse_args()
    
    input_file = args.i
    main(input_file)
    
    