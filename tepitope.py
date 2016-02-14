'''A local implementation of tepitope - an MHC class II binding predictor.'''
import os, glob

class HLAmatrix:
    '''class for loading and scoring a matrix'''

    def load_matrix(self, matrix_file):
        '''Loads and parses the matrix into memory'''
        matrix_handle = open(matrix_file)

        percentage_threshold_tmp = ''
        scoring_tmp = ''
        for line in matrix_handle:
            if len(line.split(',')[0]) <= 2:
                line = line.replace('-', '-0')
                # print line
                self.matrix_data[line.split(',')[0][0]] = map(float, (line.strip().split(',')[1:]))

            if line.split(',')[0] == 'Percent Threshold':
                percentage_threshold_tmp = line.rstrip().split(',')[1:]
            if line.split(',')[0] == 'Numerical Score':
                scoring_tmp = line.rstrip().split(',')[1:]

        for i, threshold in enumerate(percentage_threshold_tmp):
            thresh_clean = threshold.split('%')[0]
            if int(thresh_clean) == self.threshold:
                self.threshold_score = float(scoring_tmp[i])
        matrix_handle.close()

    def score(self, nonamer):
        '''Scores a nonamer based on the loaded matrix'''
        score = 0
        for i, residue in enumerate(nonamer):
            try:
                tmp_score = self.matrix_data[residue][i]
                score = score + tmp_score
            except:
                continue

        return (nonamer, score)


    def score_sequence(self, sequence):
        '''Scores an entire sequence based on the loaded matrix and threshold'''
        returned_score = []
        for i in range(len(sequence)-9):
            scored = self.score(sequence[i:i+9])
            if scored[1] >= self.threshold_score:
                returned_score.append(scored)

        return sorted(returned_score, key=lambda score: score[1], reverse=True)


    def __init__(self, csv_file, threshold):
        self.matrix_data = {}
        self.threshold = threshold
        self.threshold_score = 0
        self.matrix_file = csv_file  ## Set the matrix file variable
        self.allele = os.path.splitext(os.path.basename(csv_file))[0]    ## Set the allele name
        self.load_matrix(csv_file)   ## Load the matrix


##### Example 1 #####

sequence = "PSPPGNLRVTDVTSTSVTLSWEPPPGPITGYRVEYREAGGEWKEVTVPGSETSYTVTGLKPGTEYEFRVRAVNGAGEGPPSSVSVTT"
# print 'Sequence: %s' % (sequence)

objects = []

files = glob.glob('matrices/*.csv')
for file in files:
    objects.append(HLAmatrix(file, 5))

for matrix in objects:
    print matrix.allele
    print matrix.score_sequence(sequence)


##### Example 2 #####

# sequence = "PSPPGNLRVTDVTSTSVTLSWEPPPGPITGYRVEYREAGGEWKEVTVPGSETSYTVTGLKPGTEYEFRVRAVNGAGEGPPSSVSVTT"
# # print 'Sequence: %s' % (sequence)
# drb1_0101 = HLAmatrix('matrices/HLA-DRB1_0101.csv', 5)
# print drb1_0101.score_sequence(sequence)
# drb1_0801 = HLAmatrix('matrices/HLA-DRB1_0801.csv', 5)
# print drb1_0801.score_sequence(sequence)

