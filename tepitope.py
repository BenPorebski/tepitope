"""A local implementation of tepitope - an MHC class II binding predictor."""


class HLAmatrix:
    """class for loading and scoring a matrix"""

    def load_matrix(self, matrix_file):
        """Loads and parses the matrix into memory"""
        with open(matrix_file) as matrix_handle:
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

    def score(self, nonamer):
        """Scores a nonamer based on the loaded matrix"""
        score = 0
        for i, residue in enumerate(nonamer):
            try:
                tmp_score = self.matrix_data[residue][i]
                score = score + tmp_score
            except:
                continue

        return (nonamer, score)


    def score_sequence(self, sequence):
        """Scores an entire sequence based on the loaded matrix and threshold"""
        returned_score = []
        for i in range(len(sequence)-9):
            scored = self.score(sequence[i:i+9])
            if scored[1] >= self.threshold_score:
                returned_score.append( (scored[0],scored[1],i,i+8) )

        return sorted(returned_score, key=lambda score: score[1], reverse=True)


    def __init__(self, csv_file, threshold):
        self.matrix_data = {}
        self.threshold = threshold
        self.threshold_score = 0
        self.matrix_file = csv_file  ## Set the matrix file variable
        self.allele = os.path.splitext(os.path.basename(csv_file))[0]    ## Set the allele name
        self.load_matrix(csv_file)   ## Load the matrix


if __name__ == '__main__':

    import os, glob, sys, argparse, re

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', "--threshold", dest='threshold', type=int, default=5)
    parser.add_argument('-c', '--csv', dest='csv_filename', default=None, help='filename for CSV output')
    parser.add_argument('--no-header', dest='csv_header', action='store_false',
                        help='omit the header row from the CSV file')
    parser.add_argument('-f', '--file', '--sequence-file', dest='filename', default=None,
                        help='file which contains the sequence')
    parser.add_argument('-i', '--starting-index', dest='base_index', default=1, type=int,
                        help='the index of the first peptide in the sequence, used to offset the begin and end indexes '
                             'of the results')
    parser.add_argument("sequence", default=None, nargs='?', help='the peptide sequence to analyze, if not using -f')
    args = parser.parse_args()

    # get sequence
    if args.sequence is not None and args.filename is not None:
        print >> sys.stderr, 'You must specify only one of -f, --file, --sequence-file, or the sequence as a command-line argument'
        exit(1)
    if args.sequence is None and args.filename is None:
        print >> sys.stderr, 'You must specify -f, --file, --sequence-file, or put the sequence as a command-line argument'
        exit(1)
    if args.sequence:
        sequence = args.sequence
    else:
        with open(args.filename) as f:
            sequence = re.sub(r'\s', '', f.read())

    # generate results
    results = set()
    for matrix_file in glob.glob('matrices/*.csv'):
        matrix = HLAmatrix(matrix_file, args.threshold)
        for epitope, score, begin, end in matrix.score_sequence(sequence):
            results.add( (score, begin + args.base_index, end + args.base_index, epitope, matrix.allele) )

    results = sorted( results, key=lambda row:row[0], reverse=True )

    # output
    if args.csv_filename is not None:
        # CSV output
        import csv
        with open(args.csv_filename, 'wb') as csv_file:
            out = csv.writer(csv_file,lineterminator=os.linesep)
            if args.csv_header:
                out.writerow(('score', 'begin', 'end', 'epitope', 'allele'))
            for result in results:
                out.writerow(result)
    else:
        # stdout
        for result in results:
            print '%.2f\t%d\t%d\t%s\t%s'%result


    ##### Example 2 #####

    # sequence = "PSPPGNLRVTDVTSTSVTLSWEPPPGPITGYRVEYREAGGEWKEVTVPGSETSYTVTGLKPGTEYEFRVRAVNGAGEGPPSSVSVTT"
    # # print 'Sequence: %s' % (sequence)
    # drb1_0101 = HLAmatrix('matrices/HLA-DRB1_0101.csv', 5)
    # print drb1_0101.score_sequence(sequence)
    # drb1_0801 = HLAmatrix('matrices/HLA-DRB1_0801.csv', 5)
    # print drb1_0801.score_sequence(sequence)
