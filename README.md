
# An elementary bound on Siegel zeroes

The code in this repository verifies a calculation for L-functions associated to real Dirichlet characters, based on a method of Pintz. The code is written for Sage 8.4, in order to take advantage of its built in zetaderiv function.

## Getting Started

The code is split into two sage files depending on the parity of the desired L-functions. In the Sage terminal, navigate to the code folder and attach the desired .sage file. Sage will verify which parity has been chosen. 
```
sage: attach('even.sage')
Siegel zeroes for primitive even characters.
sage: attach('odd.sage')
Siegel zeroes for primitive odd characters.
```

## Running the verification

In the proof, we wish to show that F, a of the variables c, x, and q, is negative for some values of c and x. This is calculated by F(c, q0, q1, x).

```
sage: F(1.04100000000000, 1.30000000000000e7, 2.40000000000000e7, 10^6.83)
-0.04035111725951765
```

### Verifying Table 1

The commands to verify Table 1 are contained in data.txt in the data folder. Simply copy these into the terminal and execute.

```
sage: F(1.01100000000000, 400000, 700000, 10^5.54)
....: F(1.01700000000000, 700000, 1000000, 10^5.73)
....: F(1.02000000000000, 1000000, 1.70000000000000e6, 10^5.88)
....: F(1.02500000000000, 1.70000000000000e6, 3.10000000000000e6, 10^6.08)
....: F(1.03000000000000, 3.10000000000000e6, 6.10000000000000e6, 10^6.30)
....: F(1.03600000000000, 6.10000000000000e6, 1.30000000000000e7, 10^6.54)
....: F(1.04100000000000, 1.30000000000000e7, 2.40000000000000e7, 10^6.83)
....: F(1.04400000000000, 2.40000000000000e7, 5.40000000000000e7, 10^7.06)
....: F(1.05100000000000, 5.40000000000000e7, 1.50000000000000e8, 10^7.35)
....: F(1.05500000000000, 1.50000000000000e8, 6.20000000000000e8, 10^7.75)
....: F(1.06000000000000, 6.20000000000000e8, 4.40000000000000e9, 10^8.29)
....: F(1.07000000000000, 4.40000000000000e9, 6.40000000000000e10, 10^9.01)
....: F(1.08000000000000, 6.40000000000000e10, 2.70000000000000e12, 10^10.00)
....: F(1.09000000000000, 2.70000000000000e12, 6.20000000000000e14, 10^11.40)
....: F(1.10100000000000, 6.20000000000000e14, 2.10000000000000e18, 10^13.43)
....: F(1.20000000000000, 2.10000000000000e18, 4.40000000000000e21, 10^15.26)
....: F(1.30000000000000, 4.40000000000000e21, 1.50000000000000e24, 10^16.64)
....: F(1.35000000000000, 1.50000000000000e24, 3.10000000000000e26, 10^17.90)
....: F(1.40000000000000, 3.10000000000000e26, 2.40000000000000e28, 10^18.92)
....: F(1.42500000000000, 2.40000000000000e28, 1.50000000000000e30, 10^19.91)
....: F(1.44500000000000, 1.50000000000000e30, 100000000000000000000000000000000, 10^20.89)
....: F(1.49500000000000, 100000000000000000000000000000000, 9.10000000000000e32, 10^21.40)
....:
-0.017806884859784093
-0.020186518906101236
-0.020297219665350197
-0.01834778148933941
-0.019638045748270716
-0.005945981316594157
-0.04035111725951765
-0.023567360551733724
-0.0234476663283901
-0.023947058341990318
-0.01485948644348209
-0.006579589923060458
-0.023344215856423545
-0.03184200950716337
-0.009507322063429358
-0.006329075750147739
-0.030165661113673525
-0.016706729389835973
-0.006571171020285038
-0.0512944146355363
-0.023670093504046252
-0.026035552659428585
```

### Generating a new table

The function cq_table(q_list) generates a table of c and x values for which F(c, q0, q1, x) is negative between the values of q given in q_list. The first half of the output is the table, formatted for the tabular environment in LaTeX. The second half is the list of Sage commands to verify the table.


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thanks to Billie Thompson for the README.md template.