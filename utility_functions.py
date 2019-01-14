'''
Module:

utility_functions.py

Author:
    Jay Pedersen, UNMC, Pathology/Microbiology Department, Feb 9, 2017

Purpose:
    Utility functions used in multiple modules -- flatten (flatten a list), etc
'''

class utility_functions():

    @staticmethod
    def get_first_letter(name): # class level function, no instance (no self argument)
        # return first letter from name argument, None if no alphabetic character
        first_letter = None
        for idx,char in enumerate(name):
            next_char = name.strip()[idx].upper()
            if next_char.isalpha(): first_letter = next_char; break
        return first_letter

    # flatten: given [[a,b,c],[d,e],[f,g,h]] return [a,b,c,d,e,f,g,h]
    @staticmethod
    def flatten(list_of_lists): # class-level function
        return [item for sublist in list_of_lists for item in sublist]

# end class utility_functions
