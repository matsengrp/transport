
if 0: # process the raw yfv data from pogo github to make .tcrs files
    # https://github.com/mptouzel/pogorelyy_et_al_2018.git
    #

    #['Clone ID', 'Clone count', 'Clone fraction', 'Clonal sequence(s)', 'Clonal sequence quality(s)', 'All V hits', 'All D hits', 'All J hits', 'All C hits', 'All V alignments', 'All D alignments', 'All J alignments', 'All C alignments', 'N. Seq. FR1', 'Min. qual. FR1', 'N. Seq. CDR1', 'Min. qual. CDR1', 'N. Seq. FR2', 'Min. qual. FR2', 'N. Seq. CDR2', 'Min. qual. CDR2', 'N. Seq. FR3', 'Min. qual. FR3', 'N. Seq. CDR3', 'Min. qual. CDR3', 'N. Seq. FR4', 'Min. qual. FR4', 'AA. Seq. FR1', 'AA. Seq. CDR1', 'AA. Seq. FR2', 'AA. Seq. CDR2', 'AA. Seq. FR3', 'AA. Seq. CDR3', 'AA. Seq. FR4', 'Ref. points']

    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', \
                   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


    # make tcrs files
    files = glob('[PQS][12]*F[12]*.txt')

    topn = 1000

    for file in files:
        print(file)
        header = []
        tcrs = []
        for line in open(file,'rU'):
            l = line[:-1].split('\t')
            if not header:
                header = l[:]
                print( header )
            else:
                cdr3 = l[ header.index('AA. Seq. CDR3') ]
                badseq = len(cdr3)<=5
                for aa in cdr3:
                    if aa not in amino_acids:
                        badseq = True
                        print('bad:',cdr3)
                        break
                if badseq:
                    continue
                v_hits = l[ header.index('All V hits') ]
                vg = v_hits.split(',')[0]
                vg = vg[:vg.index('*')]+'*01'
                tcrs.append( ( vg+','+cdr3 ) )
                if len(tcrs)>=topn:
                    break

        assert len(tcrs) == topn
        outfile = '{}.top{}.tcrs'.format(file,topn)
        print('making:',outfile)
        out = open(outfile,'w')
        out.write('\n'.join(tcrs)+'\n' )
        out.close()

    exit()

