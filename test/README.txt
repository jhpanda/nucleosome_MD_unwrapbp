
usage: unwrap [parameter file]

example of [parameter file]: unwraptemp.in

what's in unwraptemp.in?
1. total numer of atoms of the nucleosome, total number of DNA atoms considered in unwrap calculations (39bp+39bp at both ends), total number of histone core heavy atoms (to calculate minimal distance from histone core, extract from "*ndx" of histone core)

#Total Nr. of atoms, DNA, base pairs & number of histone core CAs
25149 3192 6674

2. DNA heavy atom index, two columns in this section, the first one is DNA base pair number (starting from 0 to 146), the second column is index number in the nucleosome system. 
    第一列BP序列i（DNA总共有0-146），第二列是第i个BP包含的重原子index

    Note that, for the DNA base pair number, it starts from 0, ends with 38 at one end. At the other end, it starts from 77 and ends with 39. 

    Base pair number cartoon representation:
(0)------------------(38)_______________________(39)--------------------(77)

    可将DNA *itp里的原子序号和质量提取出来，如massnew.dat，执行python bp_new.py > dna.txt 得到


3. histone core heavy atom index.  
    Core index可用 python getcorendx_new.py得到，注意更改脚本里的文件名称
    可用python extract.py > core.txt得到

