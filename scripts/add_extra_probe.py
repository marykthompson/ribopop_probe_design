'''
add_extra_probe.py
Given a second csv file, add one probe to the probe set in the first csv file,
checking for heterodimer clashes.
'''
import sys
import pandas as pd
import choose_probes

def main(arglist):
    probe_file, candidate_file, outname, min_dimer_dG = sys.argv[1:]
    min_dimer_dG = float(min_dimer_dG)

    df = pd.read_csv(probe_file)
    cols = df.columns.values
    df.set_index('unique_id', inplace = True)
    df2 = pd.read_csv(candidate_file, index_col = 'unique_id')

    for i in range(0, len(df2)):
        test_df = df.append(df2.iloc[i])
        dimer_df = choose_probes.calc_dimer(test_df)
        if dimer_df.iloc[:, 0].min() > min_dimer_dG:
            test_df[['dimer_dG', 'dimer_partner']] = dimer_df
            chosen_df = test_df.copy()
            break

    chosen_df['probe_num'] = chosen_df.reset_index(drop = True).index + 1
    chosen_df.dropna(axis = 1, inplace = True)
    chosen_df.reset_index().round(2)[cols].to_csv('%s.csv' % outname, index = False)

if __name__ == '__main__':
    main(sys.argv[1:])
