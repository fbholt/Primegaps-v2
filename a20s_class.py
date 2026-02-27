import numpy as np

# Dictionaries for modeling, estimating, and enumerating constellations of lengths J=3 and J=5
J3dict = {'242': [2,4,2], '424': [4,2,4], '246': [2,4,6], '426': [4,2,6], '264': [2,6,4],
          '626': [6,2,6], '662': [6,6,2], '248': [2,4,8], '2,10,2': [2,10,2], '646': [6,4,6], '24,12': [2,4,12],
          '468': [4,6,8], '486': [4,8,6], '2,10,6': [2,10,6], '666': [6,6,6], '686': [6,8,6],
          '12,12,12': [12,12,12], '12,18,12': [12,18,12], '18,18,18': [18,18,18], 
          '30,30,30': [30,30,30], '42,42,42': [42,42,42]}

J5dict = {'42424': [4,2,4,2,4], '2,10,2,10,2': [2,10,2,10,2], '2468,10': [2,4,6,8,10],
          '66266': [6,6,2,6,6], '6,12,2,6,12': [6,12,2,6,12], '64,14,42': [6,4,14,4,2],
          '6,12,10,66': [6,12,10,6,6], '10,2,10,2,10': [10,2,10,2,10], '6,14,10,66': [6,14,10,6,6],
          '6,14,10,68': [6,14,10,6,8], '10,20,10,20,10': [10,20,10,20,10], '30,30,30,30,30': [30,30,30,30,30]}

# key arrays for consistent indexing of the keys, e.g. for legends and for DH 2d arrays
J3keyarr = ['242', '426', '24,12', '626', '662', '248', '486', '424', '264', '246', '2,10,6', '666',
            '2,10,2', '468', '2,10,6', '12,12,12', '18,18,18', '686', '42,42,42', '12,18,12', '30,30,30']

J5keyarr = ['42424', '2,10,2,10,2', '64,14,42', '66266', '6,12,2,6,12', '6,12,10,66', 
            '2468,10', '10,2,10,2,10', '6,14,10,66', '6,14,10,68', '10,20,10,20,10', '30,30,30,30,30']

# Class for admissible constellations
class AdmS:
    def __init__(self, label, gaparr):
        self.label = label
        self.gaparr = np.array(gaparr, dtype='int')

    def J(self):
        return( len(self.gaparr))

    def span(self):
        return( np.sum(self.gaparr[:]))

    def drivers(self, inputarr, cycleflag):  # returns an array of populations of driving terms for s in the input array
        i=0
        inputlen = len(inputarr)
        if (cycleflag):
            maxi = inputlen - 1
        else:
            maxi = inputlen - self.J()
        # maxi is the last possible index for a starting point for an instance of s in inputarr

        outputarr = np.zeros(30)  # length 30 would cover driving terms for spans up to 120

        while (i <= maxi):
            # does this index start a driving term for s?
            j = 0   # index in s for pattern matching
            i1 = i  # index for accumulating sums
            possible_match = True   
            
            while (j < int(self.J()) and possible_match):
                # for each component of s...
                target_sum = self.gaparr[j]
                running_sum = inputarr[i1]
                while (running_sum < target_sum and possible_match):
                    # increment the index i1....
                    if (i1 < (inputlen-1)):
                        i1 += 1
                        running_sum += inputarr[i1]
                    elif (cycleflag):
                        i1 = 0
                        running_sum += inputarr[i1]
                    else:
                        possible_match = False
                        
                if (running_sum == target_sum):
                    # this running sum matches this component of s
                    j += 1  # next element in s
                    # increment the index i1....
                    if (i1 < (inputlen-1)):
                        i1 += 1
                    elif (cycleflag):
                        i1 = 0
                    else:
                        possible_match = False
                else:
                    possible_match = False
                    
            if (possible_match):
                if (i1 > i):
                    jdex = i1-i-1
                else:
                    jdex = i1+inputlen -i -2
                if (jdex < 30) :
                    outputarr[jdex] += 1

            i += 1
        
        return outputarr

    def driversbin(self, inputarr, cycleflag):  # returns an array of populations of driving terms for s in the inputarr
        inputlen = len(inputarr)
        maxi = int((inputlen - 1) if cycleflag else (inputlen - int(self.span())-1))
        # maxi is the last possible index for a starting point for an instance of s in inputarr

        outputarr = np.zeros(30, dtype=int)  # length 30 would cover driving terms for spans up to 120

        i=0
        
        while (i <= maxi):
            # advance i to next True value
            while ((inputarr[i] != True) and (i <= maxi)):
                i += 1

            # does this index start a driving term for s?
            j = 0   # index in s for pattern matching
            possible_match = True   
            i1 = i
            
            while (j < int(self.J()) and possible_match):
                # for each component of s...
                i1 = int(i1 + self.gaparr[j]/2 ) # boolean array only records odd entries
                if (i1 >= inputlen):
                    if (cycleflag): 
                        i1 = i1 - inputlen  # loop back to beginning of cycle
                        if (not inputarr[i]): possible_match = False
                    else:
                        possible_match = False

                if (not inputarr[i1]): possible_match = False

                # print(f"i {i} i1 {i1} j {j} match {possible_match}")

                j+= 1
                    
            if (possible_match):
                if (i1 > i):
                    jdex = int(sum(inputarr[i:i1])-1)
                else:
                    jdex = int(sum(inputarr[i:])+sum(inputarr[0:i1])-1)

                # print(f"j {jdex} i {i} i1 {i1}")

                if (jdex < 30) :
                    outputarr[jdex] += 1

            if ((i % 10000000) <= 4):
                print(f"\r {i:,} of {len(inputarr):,}", end="")

            i += 1
        
        return outputarr
