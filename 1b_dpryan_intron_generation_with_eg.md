# Intronic bin generation on flattened GTF

![image](https://user-images.githubusercontent.com/68455070/127979776-c5e75e9e-ffa4-4457-a341-021fc198a609.png)

##### Table of Contents  
[Step 1: Read in flattened GTF file and view metadata](#headers)  
[Emphasis](#emphasis)  
...snip...    
<a name="headers"/>

## Step 1: Read in flattened GTF file and view metadata
```r
> gtf_test <- import.gff2("gencode.vM27.primary_assembly.annotation.dexseq.chr3.gtf")

> head(gtf_test)

GRanges object with 6 ranges and 7 metadata columns:
      seqnames          ranges strand |                       source           type     score     phase              gene_id
         <Rle>       <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>          <character>
  [1]     chr3 3069070-3069169      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA> ENSMUSG00002075012.1
  [2]     chr3 3069070-3069169      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00002075012.1
  [3]     chr3 3092718-3092817      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA> ENSMUSG00002075999.1
  [4]     chr3 3092718-3092817      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00002075999.1
  [5]     chr3 3253093-3253219      + | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA> ENSMUSG00000103416.2
  [6]     chr3 3253093-3253219      + | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000103416.2
               transcripts exonic_part_number
               <character>        <character>
  [1]                 <NA>               <NA>
  [2] ENSMUST00020182553.1                001
  [3]                 <NA>               <NA>
  [4] ENSMUST00020181724.1                001
  [5]                 <NA>               <NA>
  [6] ENSMUST00000191986.2                001
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
Get metadata.
```r
> meta_test <- elementMetadata(gtf_test)

> meta_test

DataFrame with 24904 rows and 7 columns
                            source           type     score     phase               gene_id          transcripts
                          <factor>       <factor> <numeric> <integer>           <character>          <character>
1     dexseq_prepare_annotation.py aggregate_gene        NA        NA  ENSMUSG00002075012.1                   NA
2     dexseq_prepare_annotation.py    exonic_part        NA        NA  ENSMUSG00002075012.1 ENSMUST00020182553.1
3     dexseq_prepare_annotation.py aggregate_gene        NA        NA  ENSMUSG00002075999.1                   NA
4     dexseq_prepare_annotation.py    exonic_part        NA        NA  ENSMUSG00002075999.1 ENSMUST00020181724.1
5     dexseq_prepare_annotation.py aggregate_gene        NA        NA  ENSMUSG00000103416.2                   NA
...                            ...            ...       ...       ...                   ...                  ...
24900 dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00000028173.11 ENSMUST00000198878.2
24901 dexseq_prepare_annotation.py aggregate_gene        NA        NA  ENSMUSG00000105990.2                   NA
24902 dexseq_prepare_annotation.py    exonic_part        NA        NA  ENSMUSG00000105990.2 ENSMUST00000199465.2
24903 dexseq_prepare_annotation.py aggregate_gene        NA        NA  ENSMUSG00000075903.3                   NA
24904 dexseq_prepare_annotation.py    exonic_part        NA        NA  ENSMUSG00000075903.3 ENSMUST00000101953.3
      exonic_part_number
             <character>
1                     NA
2                    001
3                     NA
4                    001
5                     NA
...                  ...
24900                025
24901                 NA
24902                001
24903                 NA
24904                001
```

Add in `intronic_part` as a factor in `gtf_test`.
```r
> elementMetadata(gtf_test)$type <- factor(elementMetadata(gtf_test)$type, levels=c(levels(elementMetadata(gtf_test)$type), "intronic_part"))

> head(elementMetadata(gtf_test)$type)

[1] aggregate_gene exonic_part    aggregate_gene exonic_part    aggregate_gene exonic_part   
Levels: aggregate_gene exonic_part intronic_part
```

Take note of rows that have NA in `exonic_part_number` column.
```r
> USE_test <- which(!is.na(elementMetadata(gtf_test)$exonic_part_number))

> USE_test

   [1]    2    4    6    8   10   11   12   13   14   15   16   17   18   19   20   21   23   25   27   29   31   33   35   37
  [25]   39   41   43   45   47   49   51   52   53   54   55   56   57   58   59   60   62   64   65   67   68   69   70   71
  [49]   73   74   75   77   78   79   80   81   82   83   84   85   86   87   88   89   90   91   92   93   94   95   96   97
  [73]   99  100  102  103  104  105  106  107  108  109  110  111  112  113  114  115  117  119  120  122  124  126  128  130
  [97]  132  134  136  138  140  141  142  143  145  147  149  151  152  153  154  155  156  157  159  161  162  164  166  167
 [121]  168  169  170  171  172  173  174  175  176  177  178  180  181  182  183  184  185  186  187  188  189  190  191  192
 [145]  193  194  196  197  198  200  202  203  204  205  206  207  209  211  213  215  217  218  219  220  221  222  223  224
 [169]  225  226  227  228  230  231  233  234  235  237  238  239  240  241  242  243  244  245  246  247  248  249  251  253
 [193]  255  256  257  258  259  261  263  264  265  267  269  271  272  273  274  275  276  277  278  279  280  281  282  283
 [217]  284  285  286  287  288  289  290  291  292  293  294  295  296  297  298  299  301  303  304  305  307  309  310  311
 [241]  312  313  314  316  318  319  320  321  322  323  324  325  327  328  329  330  331  332  333  334  335  336  337  338
 [265]  339  340  342  344  346  347  348  349  350  351  352  353  354  355  356  357  358  359  360  361  363  364  366  368
 [289]  369  371  373  374  375  376  377  378  379  380  382  383  385  386  388  390  391  392  393  394  395  396  397  398
 [313]  399  400  401  402  404  405  406  407  408  409  410  411  412  414  416  418  420  421  422  423  424  425  426  427
 [337]  428  429  430  431  433  435  437  438  439  440  441  442  443  444  445  446  447  448  449  450  451  452  453  454
 [361]  455  456  457  458  460  461  462  463  464  465  466  467  468  469  470  471  472  473  474  475  476  477  478  479
 [385]  480  481  482  483  484  485  487  488  489  490  491  492  493  494  495  496  498  499  500  501  502  503  504  505
 [409]  506  507  508  509  511  513  515  517  519  521  523  525  527  528  529  530  531  533  535  537  538  539  541  543
 [433]  545  547  549  550  552  553  554  555  556  557  558  559  560  561  562  563  564  565  566  567  568  569  570  571
 [457]  572  573  575  577  579  581  583  585  587  589  590  591  592  594  596  597  599  600  601  602  603  604  605  606
 [481]  607  609  611  612  613  614  615  616  617  618  619  620  621  622  623  624  625  626  627  628  629  630  631  632
 [505]  633  634  635  636  637  638  639  640  641  642  643  644  645  646  647  649  650  651  652  653  654  655  656  657
 [529]  658  659  661  662  663  664  665  666  667  668  669  670  671  672  673  675  677  678  679  680  681  682  683  685
 [553]  686  687  688  690  691  692  693  694  695  696  697  698  699  700  701  703  705  706  707  708  709  710  711  712
 [577]  713  714  715  717  718  719  720  721  722  723  724  725  726  727  728  730  731  732  733  735  737  738  739  740
 [601]  741  743  744  745  747  748  749  750  751  752  753  754  755  756  757  758  760  761  762  763  764  765  766  767
 [625]  768  769  770  771  773  774  776  777  778  779  780  781  782  783  784  785  786  787  788  789  790  791  792  793
 [649]  794  796  798  800  801  802  803  804  805  806  808  810  812  814  815  816  817  818  819  820  821  822  823  824
 [673]  825  827  829  831  833  835  836  837  839  841  843  844  845  846  847  848  849  850  851  852  853  854  855  856
 [697]  857  858  859  860  862  864  866  868  870  872  874  875  876  877  878  879  881  883  885  886  887  888  889  890
 [721]  892  893  894  895  897  899  901  903  905  907  908  909  910  911  912  913  915  917  918  919  920  921  922  923
 [745]  924  925  926  927  928  929  930  931  932  933  934  935  936  937  938  939  940  941  942  944  945  946  947  948
 [769]  949  950  951  952  953  954  955  956  957  958  959  960  961  962  963  964  965  966  967  968  969  970  971  972
 [793]  973  974  975  976  978  979  981  983  984  985  986  987  988  989  990  991  992  993  994  995  997  999 1000 1002
 [817] 1003 1005 1006 1007 1008 1009 1010 1011 1012 1013 1014 1015 1016 1018 1020 1021 1023 1025 1027 1028 1029 1031 1032 1034
 [841] 1036 1037 1038 1039 1040 1041 1042 1043 1044 1045 1046 1047 1048 1049 1050 1051 1052 1053 1054 1055 1056 1057 1058 1059
 [865] 1060 1061 1062 1063 1064 1065 1066 1067 1068 1069 1070 1071 1072 1073 1074 1075 1076 1077 1078 1079 1080 1081 1082 1083
 [889] 1084 1085 1086 1087 1089 1090 1091 1092 1093 1094 1095 1096 1097 1098 1099 1100 1101 1102 1103 1104 1105 1106 1107 1108
 [913] 1109 1110 1111 1112 1113 1114 1115 1116 1117 1118 1119 1120 1121 1122 1123 1124 1125 1126 1128 1129 1130 1131 1132 1133
 [937] 1134 1135 1136 1137 1138 1139 1140 1141 1142 1143 1144 1145 1146 1147 1148 1149 1150 1151 1152 1153 1154 1155 1156 1157
 [961] 1158 1159 1160 1161 1162 1163 1165 1166 1167 1168 1169 1170 1171 1172 1173 1174 1175 1176 1177 1178 1180 1181 1182 1183
 [985] 1184 1185 1186 1187 1188 1189 1190 1191 1193 1195 1196 1197 1198 1199 1200 1201
 [ reached getOption("max.print") -- omitted 21110 entries ]
```
View the `exonic_part_number` column.
```r
> elementMetadata(gtf_test)$exonic_part_number

   [1] NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011"
  [21] "012" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA   
  [41] "001" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010"
  [61] NA    "001" NA    "001" "002" NA    "001" "002" "003" "004" "005" NA    "001" "002" "003" NA    "001" "002" "003" "004"
  [81] "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" NA    "001" "002"
 [101] NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" NA    "001" NA    "001" "002"
 [121] NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001"
 [141] "002" "003" "004" NA    "001" NA    "001" NA    "001" NA    "001" "002" "003" "004" "005" "006" "007" NA    "001" NA   
 [161] "001" "002" NA    "001" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" NA    "001"
 [181] "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" NA    "001" "002" "003" NA    "001"
 [201] NA    "001" "002" "003" "004" "005" "006" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" "002" "003" "004"
 [221] "005" "006" "007" "008" "009" "010" "011" "012" NA    "001" "002" NA    "001" "002" "003" NA    "001" "002" "003" "004"
 [241] "005" "006" "007" "008" "009" "010" "011" "012" "013" NA    "001" NA    "001" NA    "001" "002" "003" "004" "005" NA   
 [261] "001" NA    "001" "002" "003" NA    "001" NA    "001" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010"
 [281] "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" "023" "024" "025" "026" "027" "028" "029" NA   
 [301] "001" NA    "001" "002" "003" NA    "001" NA    "001" "002" "003" "004" "005" "006" NA    "001" NA    "001" "002" "003"
 [321] "004" "005" "006" "007" "008" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014"
 [341] NA    "001" NA    "001" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015"
 [361] "016" NA    "001" "002" NA    "001" NA    "001" "002" NA    "001" NA    "001" "002" "003" "004" "005" "006" "007" "008"
 [381] NA    "001" "002" NA    "001" "002" NA    "001" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011"
 [401] "012" "013" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" NA    "001" NA    "001" NA    "001" NA    "001"
 [421] "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" NA    "001" NA    "001" NA    "001" "002" "003" "004"
 [441] "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" NA    "001"
 [461] "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021"
 [481] "022" "023" "024" "025" "026" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" NA    "001" "002" "003"
 [501] "004" "005" "006" "007" "008" "009" "010" "011" "012" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA   
 [521] "001" NA    "001" NA    "001" NA    "001" "002" "003" "004" "005" NA    "001" NA    "001" NA    "001" "002" "003" NA   
 [541] "001" NA    "001" NA    "001" NA    "001" NA    "001" "002" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009"
 [561] "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" NA    "001" NA    "001" NA    "001" NA   
 [581] "001" NA    "001" NA    "001" NA    "001" NA    "001" "002" "003" "004" NA    "001" NA    "001" "002" NA    "001" "002"
 [601] "003" "004" "005" "006" "007" "008" "009" NA    "001" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010"
 [621] "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" "023" "024" "025" "026" "027" "028" "029" "030"
 [641] "031" "032" "033" "034" "035" "036" "037" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" NA   
 [661] "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" NA    "001" NA    "001" "002" "003" "004"
 [681] "005" "006" "007" NA    "001" "002" "003" "004" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011"
 [701] "012" NA    "001" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" NA    "001" "002" "003" "004"
 [721] "005" "006" "007" "008" "009" "010" "011" "012" NA    "001" "002" "003" "004" NA    "001" NA    "001" "002" "003" "004"
 [741] "005" NA    "001" "002" "003" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" NA    "001"
 [761] "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" NA    "001" "002" NA    "001" "002" "003" "004" "005"
 [781] "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" NA    "001" NA    "001" NA    "001"
 [801] "002" "003" "004" "005" "006" "007" NA    "001" NA    "001" NA    "001" NA    "001" "002" "003" "004" "005" "006" "007"
 [821] "008" "009" "010" "011" "012" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" "002" "003" NA    "001" NA   
 [841] "001" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018"
 [861] NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" NA    "001" "002" "003" "004" "005" "006" NA   
 [881] "001" NA    "001" NA    "001" "002" "003" "004" "005" "006" NA    "001" "002" "003" "004" NA    "001" NA    "001" NA   
 [901] "001" NA    "001" NA    "001" NA    "001" "002" "003" "004" "005" "006" "007" NA    "001" NA    "001" "002" "003" "004"
 [921] "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" "023" "024"
 [941] "025" "026" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017"
 [961] "018" "019" "020" "021" "022" "023" "024" "025" "026" "027" "028" "029" "030" "031" "032" "033" NA    "001" "002" NA   
 [981] "001" NA    "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" NA    "001" NA    "001" "002"
 [ reached getOption("max.print") -- omitted 23904 entries ]
```
Format `exonic_part_number` and filter out NA values. View the filtered `exonic_part_number` column.

```r
> exonic_parts_test <- sprintf("%03s", elementMetadata(gtf_test)$exonic_part_number[USE_test])

> exonic_parts_test

   [1] "001" "001" "001" "001" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "001" "001" "001" "001"
  [21] "001" "001" "001" "001" "001" "001" "001" "001" "001" "001" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010"
  [41] "001" "001" "002" "001" "002" "003" "004" "005" "001" "002" "003" "001" "002" "003" "004" "005" "006" "007" "008" "009"
  [61] "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "001" "002" "001" "002" "003" "004" "005" "006"
  [81] "007" "008" "009" "010" "011" "012" "013" "014" "001" "001" "002" "001" "001" "001" "001" "001" "001" "001" "001" "001"
 [101] "001" "002" "003" "004" "001" "001" "001" "001" "002" "003" "004" "005" "006" "007" "001" "001" "002" "001" "001" "002"
 [121] "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "001" "002" "003" "004" "005" "006" "007" "008" "009"
 [141] "010" "011" "012" "013" "014" "015" "001" "002" "003" "001" "001" "002" "003" "004" "005" "006" "001" "001" "001" "001"
 [161] "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "001" "002" "001" "002" "003" "001" "002" "003"
 [181] "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "001" "001" "001" "002" "003" "004" "005" "001" "001" "002"
 [201] "003" "001" "001" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017"
 [221] "018" "019" "020" "021" "022" "023" "024" "025" "026" "027" "028" "029" "001" "001" "002" "003" "001" "001" "002" "003"
 [241] "004" "005" "006" "001" "001" "002" "003" "004" "005" "006" "007" "008" "001" "002" "003" "004" "005" "006" "007" "008"
 [261] "009" "010" "011" "012" "013" "014" "001" "001" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012"
 [281] "013" "014" "015" "016" "001" "002" "001" "001" "002" "001" "001" "002" "003" "004" "005" "006" "007" "008" "001" "002"
 [301] "001" "002" "001" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "001" "002" "003" "004"
 [321] "005" "006" "007" "008" "009" "001" "001" "001" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012"
 [341] "001" "001" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018"
 [361] "019" "020" "021" "022" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016"
 [381] "017" "018" "019" "020" "021" "022" "023" "024" "025" "026" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010"
 [401] "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "001" "001" "001" "001" "001" "001" "001" "001"
 [421] "001" "002" "003" "004" "005" "001" "001" "001" "002" "003" "001" "001" "001" "001" "001" "002" "001" "002" "003" "004"
 [441] "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" "001" "001"
 [461] "001" "001" "001" "001" "001" "001" "002" "003" "004" "001" "001" "002" "001" "002" "003" "004" "005" "006" "007" "008"
 [481] "009" "001" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018"
 [501] "019" "020" "021" "022" "023" "024" "025" "026" "027" "028" "029" "030" "031" "032" "033" "034" "035" "036" "037" "001"
 [521] "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010"
 [541] "011" "012" "013" "001" "001" "002" "003" "004" "005" "006" "007" "001" "002" "003" "004" "001" "002" "003" "004" "005"
 [561] "006" "007" "008" "009" "010" "011" "012" "001" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "001"
 [581] "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "001" "002" "003" "004" "001" "001" "002" "003" "004"
 [601] "005" "001" "002" "003" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "001" "002" "003" "004"
 [621] "005" "006" "007" "008" "009" "010" "011" "012" "001" "002" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010"
 [641] "011" "012" "013" "014" "015" "016" "017" "018" "019" "001" "001" "001" "002" "003" "004" "005" "006" "007" "001" "001"
 [661] "001" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "001" "001" "001" "001" "001" "002" "003"
 [681] "001" "001" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018"
 [701] "001" "001" "001" "001" "001" "001" "001" "002" "003" "004" "005" "006" "001" "001" "001" "002" "003" "004" "005" "006"
 [721] "001" "002" "003" "004" "001" "001" "001" "001" "001" "001" "002" "003" "004" "005" "006" "007" "001" "001" "002" "003"
 [741] "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" "023"
 [761] "024" "025" "026" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017"
 [781] "018" "019" "020" "021" "022" "023" "024" "025" "026" "027" "028" "029" "030" "031" "032" "033" "001" "002" "001" "001"
 [801] "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "001" "001" "002" "001" "002" "001" "002" "003"
 [821] "004" "005" "006" "007" "008" "009" "010" "011" "012" "001" "001" "002" "001" "001" "001" "002" "003" "001" "002" "001"
 [841] "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020"
 [861] "021" "022" "023" "024" "025" "026" "027" "028" "029" "030" "031" "032" "033" "034" "035" "036" "037" "038" "039" "040"
 [881] "041" "042" "043" "044" "045" "046" "047" "048" "049" "050" "051" "052" "001" "002" "003" "004" "005" "006" "007" "008"
 [901] "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" "023" "024" "025" "026" "027" "028"
 [921] "029" "030" "031" "032" "033" "034" "035" "036" "037" "038" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010"
 [941] "011" "012" "013" "014" "015" "016" "017" "018" "019" "020" "021" "022" "023" "024" "025" "026" "027" "028" "029" "030"
 [961] "031" "032" "033" "034" "035" "036" "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "013" "014"
 [981] "001" "002" "003" "004" "005" "006" "007" "008" "009" "010" "011" "012" "001" "001" "002" "003" "004" "005" "006" "007"
 [ reached getOption("max.print") -- omitted 21110 entries ]
 ```
Change the `exonic_part_number` column to character data type.

```r
> elementMetadata(gtf_test)$exonic_part_number <- as.character(elementMetadata(gtf_test)$exonic_part_number)

> head(elementMetadata(gtf_test))

DataFrame with 6 rows and 7 columns
                        source           type     score     phase              gene_id          transcripts exonic_part_number
                      <factor>       <factor> <numeric> <integer>          <character>          <character>        <character>
1 dexseq_prepare_annotation.py aggregate_gene        NA        NA ENSMUSG00002075012.1                   NA                 NA
2 dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00002075012.1 ENSMUST00020182553.1                001
3 dexseq_prepare_annotation.py aggregate_gene        NA        NA ENSMUSG00002075999.1                   NA                 NA
4 dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00002075999.1 ENSMUST00020181724.1                001
5 dexseq_prepare_annotation.py aggregate_gene        NA        NA ENSMUSG00000103416.2                   NA                 NA
6 dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00000103416.2 ENSMUST00000191986.2                001
```
Replace `exonic_part_number` column with filtered version, without NAs.
```r
> elementMetadata(gtf_test)$exonic_part_number[USE_test] <- exonic_parts_test

> elementMetadata(gtf_test)

DataFrame with 24904 rows and 7 columns
                            source           type     score     phase               gene_id          transcripts
                          <factor>       <factor> <numeric> <integer>           <character>          <character>
1     dexseq_prepare_annotation.py aggregate_gene        NA        NA  ENSMUSG00002075012.1                   NA
2     dexseq_prepare_annotation.py    exonic_part        NA        NA  ENSMUSG00002075012.1 ENSMUST00020182553.1
3     dexseq_prepare_annotation.py aggregate_gene        NA        NA  ENSMUSG00002075999.1                   NA
4     dexseq_prepare_annotation.py    exonic_part        NA        NA  ENSMUSG00002075999.1 ENSMUST00020181724.1
5     dexseq_prepare_annotation.py aggregate_gene        NA        NA  ENSMUSG00000103416.2                   NA
...                            ...            ...       ...       ...                   ...                  ...
24900 dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00000028173.11 ENSMUST00000198878.2
24901 dexseq_prepare_annotation.py aggregate_gene        NA        NA  ENSMUSG00000105990.2                   NA
24902 dexseq_prepare_annotation.py    exonic_part        NA        NA  ENSMUSG00000105990.2 ENSMUST00000199465.2
24903 dexseq_prepare_annotation.py aggregate_gene        NA        NA  ENSMUSG00000075903.3                   NA
24904 dexseq_prepare_annotation.py    exonic_part        NA        NA  ENSMUSG00000075903.3 ENSMUST00000101953.3
      exonic_part_number
             <character>
1                     NA
2                    001
3                     NA
4                    001
5                     NA
...                  ...
24900                025
24901                 NA
24902                001
24903                 NA
24904                001

> gtf_test

GRanges object with 24904 ranges and 7 metadata columns:
          seqnames              ranges strand |                       source           type     score     phase
             <Rle>           <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>
      [1]     chr3     3069070-3069169      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA>
      [2]     chr3     3069070-3069169      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
      [3]     chr3     3092718-3092817      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA>
      [4]     chr3     3092718-3092817      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
      [5]     chr3     3253093-3253219      + | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA>
      ...      ...                 ...    ... .                          ...            ...       ...       ...
  [24900]     chr3 159644091-159644300      + | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
  [24901]     chr3 159581483-159583853      + | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA>
  [24902]     chr3 159581483-159583853      + | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
  [24903]     chr3 159621892-159621998      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA>
  [24904]     chr3 159621892-159621998      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
                        gene_id          transcripts exonic_part_number
                    <character>          <character>        <character>
      [1]  ENSMUSG00002075012.1                 <NA>               <NA>
      [2]  ENSMUSG00002075012.1 ENSMUST00020182553.1                001
      [3]  ENSMUSG00002075999.1                 <NA>               <NA>
      [4]  ENSMUSG00002075999.1 ENSMUST00020181724.1                001
      [5]  ENSMUSG00000103416.2                 <NA>               <NA>
      ...                   ...                  ...                ...
  [24900] ENSMUSG00000028173.11 ENSMUST00000198878.2                025
  [24901]  ENSMUSG00000105990.2                 <NA>               <NA>
  [24902]  ENSMUSG00000105990.2 ENSMUST00000199465.2                001
  [24903]  ENSMUSG00000075903.3                 <NA>               <NA>
  [24904]  ENSMUSG00000075903.3 ENSMUST00000101953.3                001
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
## Step 2: Manipulate the GTF file into a GRangesList, with each gene as a GRanges object. 

Within each gene, we have intervals of exons as GRanges ranges. 

```r
> grl_test <- split(gtf_test, elementMetadata(gtf_test)$gene_id)

> grl_test

GRangesList object of length 2794:
$ENSMUSG00000000001.5
GRanges object with 10 ranges and 7 metadata columns:
       seqnames              ranges strand |                       source           type     score     phase
          <Rle>           <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>
   [1]     chr3 108014596-108053462      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA>
   [2]     chr3 108014596-108016632      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [3]     chr3 108016719-108016928      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [4]     chr3 108019251-108019404      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [5]     chr3 108019789-108019918      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [6]     chr3 108023079-108023207      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [7]     chr3 108025617-108025774      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [8]     chr3 108030858-108030999      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [9]     chr3 108031111-108031153      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
  [10]     chr3 108053204-108053462      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
                    gene_id          transcripts exonic_part_number
                <character>          <character>        <character>
   [1] ENSMUSG00000000001.5                 <NA>               <NA>
   [2] ENSMUSG00000000001.5 ENSMUST00000000001.5                001
   [3] ENSMUSG00000000001.5 ENSMUST00000000001.5                002
   [4] ENSMUSG00000000001.5 ENSMUST00000000001.5                003
   [5] ENSMUSG00000000001.5 ENSMUST00000000001.5                004
   [6] ENSMUSG00000000001.5 ENSMUST00000000001.5                005
   [7] ENSMUSG00000000001.5 ENSMUST00000000001.5                006
   [8] ENSMUSG00000000001.5 ENSMUST00000000001.5                007
   [9] ENSMUSG00000000001.5 ENSMUST00000000001.5                008
  [10] ENSMUSG00000000001.5 ENSMUST00000000001.5                009
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

...
<2793 more elements>
```
## Step 3: Start to generate intronic bins within GTF file.

### Breakdown of Function 1: add_introns
The `add_introns()` function takes in a GRanges object containing only exonic regions, and inserts intronic parts between the existing exonic regions. 

Whole function:
```r
add_introns <- function(gr) {
  exons <- gr[which(elementMetadata(gr)$type=="exonic_part"),] # load in exonic parts types
  if(length(exons) > 1) {# exons=a list of exons 
    seqname <- seqnames(exons)[-1]# all but the 1st element
    starts <- end(exons)+1# end(x) grabs the last value of an interval (we are starting at the end of exons) Eg. 124321 + 1
    starts <- starts[-length(starts)] # gets the length of starts (1)
    ends <- start(exons)-1# get the start value of exon interval +1
    ends <- ends[-1]# all but the first element
    bounds <- IRanges(start=starts, end=ends) # defining integer ranges
    strand <- strand(exons)[-1] # all but the first element
    # GRanges object is a collection of genomic features 
    introns <- GRanges(seqnames=seqname, ranges=bounds, strand=strand) # constructing a new GRanges object
    intron_ids <- sprintf("%03i", c(1:length(introns))) # naming the introns (intron_ids) 
    
    # Remove 0-width introns based on their ranges
    DISCARD <- which(width(introns) <= 0) 
    # which() returns a vector of introns that contains width <= 0
    # the vector of introns to discard
    # QNS: WHY will there be introns that are negative width?
    if(length(DISCARD) > 0) { # if there are things to discard
      introns <- introns[-DISCARD] # grabs all but the DISCARD elements
      intron_ids <- intron_ids[-DISCARD] # grabs all but the DISCARD id
    }
    if(length(introns) > 0) {# if there are introns to create
      # create the meta-data
      df <- as.data.frame(elementMetadata(exons))# load in the current exons as a meta-dataframe where exons <- gr[which(elementMetadata(gr)$type=="exonic_part"),]
      nrows <- length(introns) # number of new rows to add for introns
      
      # i think concatenating df (exons) with empty rows for introns
      metadf <- df[1:nrows,] # does this need to deal with gene_id and transcripts differently?
      
      # basically just transforming the gene_id and transcripts cols to character
      metadf <- transform(metadf, gene_id=as.character(gene_id), transcripts=as.character(transcripts))
      
      # replicate "NA" nrows times (i think this is added as rows)
      metadf$transcripts <- as.character(c(rep(NA, nrows)))
      
      # replicate "intronic_part" nrows times under the $type column --> then factor the $type column in metadf
      metadf$type <- factor(c(rep("intronic_part", nrows)), levels=levels(metadf$type))
      
      # appending intron_ids to $exonic_part_number column in metadf
      metadf$exonic_part_number <- intron_ids # we created the intron_id earlier
      
      # similar to mcols function (get or set metadata columns)
      #i think i load in metadf into introns metadata
      elementMetadata(introns) <- metadf
      
      # Merge the GRanges containing exons and introns
      # introns is its own metadata
      gr <- append(gr, introns)
      gr <- gr[order(start(gr), elementMetadata(gr)$type),] # re-sort from start of gr, according to type (exonic_part or intronic_part)
    }
  }
  return(gr)
}
```

Function breakdown:

We will use the first gene, `ENSMUSG00000000001.5`, in `grl_test` as an example. 
```r
> grl_test$ENSMUSG00000000001.5

GRanges object with 10 ranges and 7 metadata columns:
       seqnames              ranges strand |                       source           type     score     phase
          <Rle>           <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>
   [1]     chr3 108014596-108053462      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA>
   [2]     chr3 108014596-108016632      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [3]     chr3 108016719-108016928      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [4]     chr3 108019251-108019404      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [5]     chr3 108019789-108019918      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [6]     chr3 108023079-108023207      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [7]     chr3 108025617-108025774      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [8]     chr3 108030858-108030999      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [9]     chr3 108031111-108031153      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
  [10]     chr3 108053204-108053462      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
                    gene_id          transcripts exonic_part_number
                <character>          <character>        <character>
   [1] ENSMUSG00000000001.5                 <NA>               <NA>
   [2] ENSMUSG00000000001.5 ENSMUST00000000001.5                001
   [3] ENSMUSG00000000001.5 ENSMUST00000000001.5                002
   [4] ENSMUSG00000000001.5 ENSMUST00000000001.5                003
   [5] ENSMUSG00000000001.5 ENSMUST00000000001.5                004
   [6] ENSMUSG00000000001.5 ENSMUST00000000001.5                005
   [7] ENSMUSG00000000001.5 ENSMUST00000000001.5                006
   [8] ENSMUSG00000000001.5 ENSMUST00000000001.5                007
   [9] ENSMUSG00000000001.5 ENSMUST00000000001.5                008
  [10] ENSMUSG00000000001.5 ENSMUST00000000001.5                009
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
  
> elementMetadata(grl_test$ENSMUSG00000000001.5)

DataFrame with 10 rows and 7 columns
                         source           type     score     phase              gene_id          transcripts exonic_part_number
                       <factor>       <factor> <numeric> <integer>          <character>          <character>        <character>
1  dexseq_prepare_annotation.py aggregate_gene        NA        NA ENSMUSG00000000001.5                   NA                 NA
2  dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00000000001.5 ENSMUST00000000001.5                001
3  dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00000000001.5 ENSMUST00000000001.5                002
4  dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00000000001.5 ENSMUST00000000001.5                003
5  dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00000000001.5 ENSMUST00000000001.5                004
6  dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00000000001.5 ENSMUST00000000001.5                005
7  dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00000000001.5 ENSMUST00000000001.5                006
8  dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00000000001.5 ENSMUST00000000001.5                007
9  dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00000000001.5 ENSMUST00000000001.5                008
10 dexseq_prepare_annotation.py    exonic_part        NA        NA ENSMUSG00000000001.5 ENSMUST00000000001.5                009

> exons <- grl_test$ENSMUSG00000000001.5[which(elementMetadata(grl_test$ENSMUSG00000000001.5)$type=="exonic_part"),]

> exons

GRanges object with 9 ranges and 7 metadata columns:
      seqnames              ranges strand |                       source        type     score     phase              gene_id
         <Rle>           <IRanges>  <Rle> |                     <factor>    <factor> <numeric> <integer>          <character>
  [1]     chr3 108014596-108016632      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5
  [2]     chr3 108016719-108016928      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5
  [3]     chr3 108019251-108019404      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5
  [4]     chr3 108019789-108019918      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5
  [5]     chr3 108023079-108023207      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5
  [6]     chr3 108025617-108025774      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5
  [7]     chr3 108030858-108030999      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5
  [8]     chr3 108031111-108031153      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5
  [9]     chr3 108053204-108053462      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5
               transcripts exonic_part_number
               <character>        <character>
  [1] ENSMUST00000000001.5                001
  [2] ENSMUST00000000001.5                002
  [3] ENSMUST00000000001.5                003
  [4] ENSMUST00000000001.5                004
  [5] ENSMUST00000000001.5                005
  [6] ENSMUST00000000001.5                006
  [7] ENSMUST00000000001.5                007
  [8] ENSMUST00000000001.5                008
  [9] ENSMUST00000000001.5                009
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
  
> length(exons)

[1] 9
```

Get all exons but the last entry (since the last exon entry does not require a preceding intronic region).
```r
> seqnames(exons)[-1]

factor-Rle of length 8 with 1 run
  Lengths:    8
  Values : chr3
Levels(1): chr3

> seqname <- seqnames(exons)[-1]
```

Generate the starts of all the intronic region (which is simply the exonic end + 1).
```r
> end(exons)

[1] 108016632 108016928 108019404 108019918 108023207 108025774 108030999 108031153 108053462

> end(exons)+1

[1] 108016633 108016929 108019405 108019919 108023208 108025775 108031000 108031154 108053463

> starts <- end(exons)+1
```

Get all the intronic starts except for the last intronic region. 
```r
> starts[-length(starts)]

[1] 108016633 108016929 108019405 108019919 108023208 108025775 108031000 108031154

> starts <- starts[-length(starts)]
```

Now we generate the ends of all the intronic regions. 
```r
> start(exons)

[1] 108014596 108016719 108019251 108019789 108023079 108025617 108030858 108031111 108053204

> start(exons)-1

[1] 108014595 108016718 108019250 108019788 108023078 108025616 108030857 108031110 108053203

> ends <- start(exons)-1
```

Get all but the first end of intronic region.
```r
> ends[-1]

[1] 108016718 108019250 108019788 108023078 108025616 108030857 108031110 108053203

> ends <- ends[-1]
```

Generate an IRanges object from the starts and ends of the intronic regions.
```r
> bounds <- IRanges(start=starts, end=ends)

> bounds

IRanges object with 8 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1] 108016633 108016718        86
  [2] 108016929 108019250      2322
  [3] 108019405 108019788       384
  [4] 108019919 108023078      3160
  [5] 108023208 108025616      2409
  [6] 108025775 108030857      5083
  [7] 108031000 108031110       111
  [8] 108031154 108053203     22050
```

Getting the strand information for all but the last exonic region.
```r
> strand(exons)[-1]

factor-Rle of length 8 with 1 run
  Lengths: 8
  Values : -
Levels(3): + - *

> strand <- strand(exons)[-1]
```

Construct a new GRanges object, which holds all the **intronic** intervals sandwiched by the exonic intervals. By comparing `introns` to `exons` GRanges objects, we can see how the exonic regions are interspersed by intronic bins. 
```r
> introns <- GRanges(seqnames=seqname, ranges=bounds, strand=strand)

> introns

GRanges object with 8 ranges and 0 metadata columns:
      seqnames              ranges strand
         <Rle>           <IRanges>  <Rle>
  [1]     chr3 108016633-108016718      -
  [2]     chr3 108016929-108019250      -
  [3]     chr3 108019405-108019788      -
  [4]     chr3 108019919-108023078      -
  [5]     chr3 108023208-108025616      -
  [6]     chr3 108025775-108030857      -
  [7]     chr3 108031000-108031110      -
  [8]     chr3 108031154-108053203      -
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
  
> exons

GRanges object with 9 ranges and 7 metadata columns:
      seqnames              ranges strand |                       source        type     score     phase              gene_id          transcripts
         <Rle>           <IRanges>  <Rle> |                     <factor>    <factor> <numeric> <integer>          <character>          <character>
  [1]     chr3 108014596-108016632      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5 ENSMUST00000000001.5
  [2]     chr3 108016719-108016928      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5 ENSMUST00000000001.5
  [3]     chr3 108019251-108019404      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5 ENSMUST00000000001.5
  [4]     chr3 108019789-108019918      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5 ENSMUST00000000001.5
  [5]     chr3 108023079-108023207      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5 ENSMUST00000000001.5
  [6]     chr3 108025617-108025774      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5 ENSMUST00000000001.5
  [7]     chr3 108030858-108030999      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5 ENSMUST00000000001.5
  [8]     chr3 108031111-108031153      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5 ENSMUST00000000001.5
  [9]     chr3 108053204-108053462      - | dexseq_prepare_annotation.py exonic_part      <NA>      <NA> ENSMUSG00000000001.5 ENSMUST00000000001.5
      exonic_part_number
             <character>
  [1]                001
  [2]                002
  [3]                003
  [4]                004
  [5]                005
  [6]                006
  [7]                007
  [8]                008
  [9]                009
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

We now want to "annotate" the intronic bins with intron IDs (arbitrarily numbered 1, 2, 3, etc.)

```r
> sprintf("%03i", c(1:length(introns)))

[1] "001" "002" "003" "004" "005" "006" "007" "008"

> intron_ids <- sprintf("%03i", c(1:length(introns)))
```

To account for a special case where the intronic bin is 0 bp large (which can occur when 2 exons are consecutively next to each other), we get rid of the 0bp intronic bins. 

```r
> which(width(introns) <= 0)

integer(0)

> width(introns) <= 0

[1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

> DISCARD <- which(width(introns) <= 0)

if(length(DISCARD) > 0) { # if there are things to discard
  introns <- introns[-DISCARD] # grabs all but the DISCARD elements and updates introns
  intron_ids <- intron_ids[-DISCARD] # grabs all but the DISCARD id and updates intron_ids
}
```

If there are intronic bins created for that gene (i.e. `if(length(introns) > 0)`), we proceed as follows:

```r
> as.data.frame(elementMetadata(exons))

                        source        type score phase              gene_id          transcripts exonic_part_number
1 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                001
2 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                002
3 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                003
4 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                004
5 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                005
6 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                006
7 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                007
8 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                008
9 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                009
> df <- as.data.frame(elementMetadata(exons))
```

Get the number of additional rows to add for intronic bins.

```r
> nrows <- length(introns)

> length(introns) 

[1] 8
```

Get the current exonic bins based on the number of intronic bins.

```r
> df[1:nrows,]

                        source        type score phase              gene_id          transcripts exonic_part_number
1 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                001
2 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                002
3 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                003
4 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                004
5 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                005
6 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                006
7 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                007
8 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5 ENSMUST00000000001.5                008

> metadf <- df[1:nrows,]
```

Transform the `gene_id` and `transcripts` columns to `character` type.

```r
> class(metadf$gene_id)

[1] "character"

> class(metadf$transcripts)

[1] "character"

> metadf <- transform(metadf, gene_id=as.character(gene_id), transcripts=as.character(transcripts))
```

Fill the `transcripts` values to all `NA`, since we are creating a dataframe for intronic bins.

```r
> as.character(c(rep(NA, nrows)))

[1] NA NA NA NA NA NA NA NA

> metadf$transcripts <- as.character(c(rep(NA, nrows)))

> metadf

                        source        type score phase              gene_id transcripts exonic_part_number
1 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5        <NA>                001
2 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5        <NA>                002
3 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5        <NA>                003
4 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5        <NA>                004
5 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5        <NA>                005
6 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5        <NA>                006
7 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5        <NA>                007
8 dexseq_prepare_annotation.py exonic_part    NA    NA ENSMUSG00000000001.5        <NA>                008
```

Fill in the `type` column for the dataframe as `intronic_part`.

```r
> c(rep("intronic_part", nrows))

[1] "intronic_part" "intronic_part" "intronic_part" "intronic_part" "intronic_part" "intronic_part" "intronic_part" "intronic_part"

> factor(c(rep("intronic_part", nrows)), levels=levels(metadf$type))

[1] intronic_part intronic_part intronic_part intronic_part intronic_part intronic_part intronic_part intronic_part

Levels: aggregate_gene exonic_part intronic_part

> metadf$type <- factor(c(rep("intronic_part", nrows)), levels=levels(metadf$type))

> metadf

                        source          type score phase              gene_id transcripts exonic_part_number
1 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                001
2 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                002
3 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                003
4 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                004
5 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                005
6 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                006
7 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                007
8 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                008
```

Replace `exonic_part_number` with the intron IDs that we've created earlier.

```r
> intron_ids

[1] "001" "002" "003" "004" "005" "006" "007" "008"

> metadf$exonic_part_number <- intron_ids

> metadf

                        source          type score phase              gene_id transcripts exonic_part_number
1 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                001
2 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                002
3 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                003
4 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                004
5 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                005
6 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                006
7 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                007
8 dexseq_prepare_annotation.py intronic_part    NA    NA ENSMUSG00000000001.5        <NA>                008
```

Now that the intron metadata has been created, load it into the `introns` GRanges we've created earlier.

```r
> introns

GRanges object with 8 ranges and 0 metadata columns:
      seqnames              ranges strand
         <Rle>           <IRanges>  <Rle>
  [1]     chr3 108016633-108016718      -
  [2]     chr3 108016929-108019250      -
  [3]     chr3 108019405-108019788      -
  [4]     chr3 108019919-108023078      -
  [5]     chr3 108023208-108025616      -
  [6]     chr3 108025775-108030857      -
  [7]     chr3 108031000-108031110      -
  [8]     chr3 108031154-108053203      -
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
  
> elementMetadata(introns) <- metadf

> introns

GRanges object with 8 ranges and 7 metadata columns:
      seqnames              ranges strand |                       source          type     score     phase              gene_id transcripts
         <Rle>           <IRanges>  <Rle> |                     <factor>      <factor> <numeric> <integer>          <character> <character>
  [1]     chr3 108016633-108016718      - | dexseq_prepare_annotation.py intronic_part      <NA>      <NA> ENSMUSG00000000001.5        <NA>
  [2]     chr3 108016929-108019250      - | dexseq_prepare_annotation.py intronic_part      <NA>      <NA> ENSMUSG00000000001.5        <NA>
  [3]     chr3 108019405-108019788      - | dexseq_prepare_annotation.py intronic_part      <NA>      <NA> ENSMUSG00000000001.5        <NA>
  [4]     chr3 108019919-108023078      - | dexseq_prepare_annotation.py intronic_part      <NA>      <NA> ENSMUSG00000000001.5        <NA>
  [5]     chr3 108023208-108025616      - | dexseq_prepare_annotation.py intronic_part      <NA>      <NA> ENSMUSG00000000001.5        <NA>
  [6]     chr3 108025775-108030857      - | dexseq_prepare_annotation.py intronic_part      <NA>      <NA> ENSMUSG00000000001.5        <NA>
  [7]     chr3 108031000-108031110      - | dexseq_prepare_annotation.py intronic_part      <NA>      <NA> ENSMUSG00000000001.5        <NA>
  [8]     chr3 108031154-108053203      - | dexseq_prepare_annotation.py intronic_part      <NA>      <NA> ENSMUSG00000000001.5        <NA>
      exonic_part_number
             <character>
  [1]                001
  [2]                002
  [3]                003
  [4]                004
  [5]                005
  [6]                006
  [7]                007
  [8]                008
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
Now the intronic bin GRanges are suitable for combining with the exons Granges. We see that the final GRanges object has both `exonic_part` (exonic bins) and `intronic_part` (intronic bins). 

```r
> grl_test$ENSMUSG00000000001.5

GRanges object with 10 ranges and 7 metadata columns:
       seqnames              ranges strand |                       source           type     score     phase              gene_id
          <Rle>           <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>          <character>
   [1]     chr3 108014596-108053462      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA> ENSMUSG00000000001.5
   [2]     chr3 108014596-108016632      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [3]     chr3 108016719-108016928      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [4]     chr3 108019251-108019404      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [5]     chr3 108019789-108019918      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [6]     chr3 108023079-108023207      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [7]     chr3 108025617-108025774      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [8]     chr3 108030858-108030999      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [9]     chr3 108031111-108031153      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
  [10]     chr3 108053204-108053462      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
                transcripts exonic_part_number
                <character>        <character>
   [1]                 <NA>               <NA>
   [2] ENSMUST00000000001.5                001
   [3] ENSMUST00000000001.5                002
   [4] ENSMUST00000000001.5                003
   [5] ENSMUST00000000001.5                004
   [6] ENSMUST00000000001.5                005
   [7] ENSMUST00000000001.5                006
   [8] ENSMUST00000000001.5                007
   [9] ENSMUST00000000001.5                008
  [10] ENSMUST00000000001.5                009
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
  
> append(grl_test$ENSMUSG00000000001.5, introns)

GRanges object with 18 ranges and 7 metadata columns:
       seqnames              ranges strand |                       source           type     score     phase              gene_id
          <Rle>           <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>          <character>
   [1]     chr3 108014596-108053462      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA> ENSMUSG00000000001.5
   [2]     chr3 108014596-108016632      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [3]     chr3 108016719-108016928      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [4]     chr3 108019251-108019404      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [5]     chr3 108019789-108019918      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   ...      ...                 ...    ... .                          ...            ...       ...       ...                  ...
  [14]     chr3 108019919-108023078      - | dexseq_prepare_annotation.py  intronic_part      <NA>      <NA> ENSMUSG00000000001.5
  [15]     chr3 108023208-108025616      - | dexseq_prepare_annotation.py  intronic_part      <NA>      <NA> ENSMUSG00000000001.5
  [16]     chr3 108025775-108030857      - | dexseq_prepare_annotation.py  intronic_part      <NA>      <NA> ENSMUSG00000000001.5
  [17]     chr3 108031000-108031110      - | dexseq_prepare_annotation.py  intronic_part      <NA>      <NA> ENSMUSG00000000001.5
  [18]     chr3 108031154-108053203      - | dexseq_prepare_annotation.py  intronic_part      <NA>      <NA> ENSMUSG00000000001.5
                transcripts exonic_part_number
                <character>        <character>
   [1]                 <NA>               <NA>
   [2] ENSMUST00000000001.5                001
   [3] ENSMUST00000000001.5                002
   [4] ENSMUST00000000001.5                003
   [5] ENSMUST00000000001.5                004
   ...                  ...                ...
  [14]                 <NA>                004
  [15]                 <NA>                005
  [16]                 <NA>                006
  [17]                 <NA>                007
  [18]                 <NA>                008
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

We sort/order the final GRanges output (containing both exonic and intronic bins) according to:

1. Their start coordinate
2. Their type (aggregate_gene --> exonic_part --> intronic_part)

```r
> elementMetadata(grl_test$ENSMUSG00000000001.5)$type

 [1] aggregate_gene exonic_part    exonic_part    exonic_part    exonic_part    exonic_part    exonic_part    exonic_part    exonic_part   
[10] exonic_part   
Levels: aggregate_gene exonic_part intronic_part

> start(grl_test$ENSMUSG00000000001.5)

 [1] 108014596 108014596 108016719 108019251 108019789 108023079 108025617 108030858 108031111 108053204
 
> grl_test$ENSMUSG00000000001.5[order(start(grl_test$ENSMUSG00000000001.5), elementMetadata(grl_test$ENSMUSG00000000001.5)$type),]

GRanges object with 10 ranges and 7 metadata columns:
       seqnames              ranges strand |                       source           type     score     phase              gene_id
          <Rle>           <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>          <character>
   [1]     chr3 108014596-108053462      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA> ENSMUSG00000000001.5
   [2]     chr3 108014596-108016632      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [3]     chr3 108016719-108016928      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [4]     chr3 108019251-108019404      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [5]     chr3 108019789-108019918      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [6]     chr3 108023079-108023207      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [7]     chr3 108025617-108025774      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [8]     chr3 108030858-108030999      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
   [9]     chr3 108031111-108031153      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
  [10]     chr3 108053204-108053462      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000000001.5
                transcripts exonic_part_number
                <character>        <character>
   [1]                 <NA>               <NA>
   [2] ENSMUST00000000001.5                001
   [3] ENSMUST00000000001.5                002
   [4] ENSMUST00000000001.5                003
   [5] ENSMUST00000000001.5                004
   [6] ENSMUST00000000001.5                005
   [7] ENSMUST00000000001.5                006
   [8] ENSMUST00000000001.5                007
   [9] ENSMUST00000000001.5                008
  [10] ENSMUST00000000001.5                009
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

The preceding `add_introns` function is applied to **all** the genes within the original GRangesList, with each gene being a GRanges object within the GRangesList. Once that is done, we can see that each gene has counting bins not only for exonic parts, but also intronic parts. There is nothing "special" about these intronic bins - they are simply coordinates of the genomic intervals between existing exonic parts. 

```r
# Original GRangesList
> grl_test

GRangesList object of length 2794:
$ENSMUSG00000000001.5
GRanges object with 10 ranges and 7 metadata columns:
       seqnames              ranges strand |                       source           type     score     phase
          <Rle>           <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>
   [1]     chr3 108014596-108053462      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA>
   [2]     chr3 108014596-108016632      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [3]     chr3 108016719-108016928      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [4]     chr3 108019251-108019404      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [5]     chr3 108019789-108019918      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [6]     chr3 108023079-108023207      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [7]     chr3 108025617-108025774      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [8]     chr3 108030858-108030999      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
   [9]     chr3 108031111-108031153      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
  [10]     chr3 108053204-108053462      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA>
                    gene_id          transcripts exonic_part_number
                <character>          <character>        <character>
   [1] ENSMUSG00000000001.5                 <NA>               <NA>
   [2] ENSMUSG00000000001.5 ENSMUST00000000001.5                001
   [3] ENSMUSG00000000001.5 ENSMUST00000000001.5                002
   [4] ENSMUSG00000000001.5 ENSMUST00000000001.5                003
   [5] ENSMUSG00000000001.5 ENSMUST00000000001.5                004
   [6] ENSMUSG00000000001.5 ENSMUST00000000001.5                005
   [7] ENSMUSG00000000001.5 ENSMUST00000000001.5                006
   [8] ENSMUSG00000000001.5 ENSMUST00000000001.5                007
   [9] ENSMUSG00000000001.5 ENSMUST00000000001.5                008
  [10] ENSMUSG00000000001.5 ENSMUST00000000001.5                009
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

...
<2793 more elements>
# Result of applying with_introns function to each gene within the GRangesList
> with_introns_grl_test

GRangesList object of length 2794:
$ENSMUSG00000000001.5
GRanges object with 18 ranges and 7 metadata columns:
       seqnames              ranges strand |                       source           type     score     phase
          <Rle>           <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>
   [1]     chr3 108014596-108053462      - | dexseq_prepare_annotation.py aggregate_gene        NA      <NA>
   [2]     chr3 108014596-108016632      - | dexseq_prepare_annotation.py exonic_part           NA      <NA>
   [3]     chr3 108016633-108016718      - | dexseq_prepare_annotation.py intronic_part         NA      <NA>
   [4]     chr3 108016719-108016928      - | dexseq_prepare_annotation.py exonic_part           NA      <NA>
   [5]     chr3 108016929-108019250      - | dexseq_prepare_annotation.py intronic_part         NA      <NA>
   ...      ...                 ...    ... .                          ...            ...       ...       ...
  [14]     chr3 108030858-108030999      - | dexseq_prepare_annotation.py  exonic_part          NA      <NA>
  [15]     chr3 108031000-108031110      - | dexseq_prepare_annotation.py  intronic_part        NA      <NA>
  [16]     chr3 108031111-108031153      - | dexseq_prepare_annotation.py  exonic_part          NA      <NA>
  [17]     chr3 108031154-108053203      - | dexseq_prepare_annotation.py  intronic_part        NA      <NA>
  [18]     chr3 108053204-108053462      - | dexseq_prepare_annotation.py  exonic_part          NA      <NA>
                    gene_id          transcripts exonic_part_number
                <character>          <character>        <character>
   [1] ENSMUSG00000000001.5                 <NA>               <NA>
   [2] ENSMUSG00000000001.5 ENSMUST00000000001.5                001
   [3] ENSMUSG00000000001.5                 <NA>                001
   [4] ENSMUSG00000000001.5 ENSMUST00000000001.5                002
   [5] ENSMUSG00000000001.5                 <NA>                002
   ...                  ...                  ...                ...
  [14] ENSMUSG00000000001.5 ENSMUST00000000001.5                007
  [15] ENSMUSG00000000001.5                 <NA>                007
  [16] ENSMUSG00000000001.5 ENSMUST00000000001.5                008
  [17] ENSMUSG00000000001.5                 <NA>                008
  [18] ENSMUSG00000000001.5 ENSMUST00000000001.5                009
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

...
<2793 more elements>
```

Before moving on, we do some simple sorting of the intervals based on chromosome and start coordinate, to facilitate downstream steps. To do this, we first grab the chromosome and start coordinates of all the intervals. We can see the difference between a sorted and unsorted output. In the unsorted output, the GRangesList is sorted by ENSMUSG gene name, while the sorted output is sorted based on actual coordinates. 

```r
> chroms <- sapply(with_introns_grl_test, function(x) as.factor(seqnames(x))[1])

> head(chroms)

 ENSMUSG00000000001.5 ENSMUSG00000000339.15 ENSMUSG00000000340.11 ENSMUSG00000000563.18 ENSMUSG00000000794.10 ENSMUSG00000001016.13 
                 chr3                  chr3                  chr3                  chr3                  chr3                  chr3 
Levels: chr3

> starts <- sapply(with_introns_grl_test, function(x) start(x)[1])

> head(starts)

 ENSMUSG00000000001.5 ENSMUSG00000000339.15 ENSMUSG00000000340.11 ENSMUSG00000000563.18 ENSMUSG00000000794.10 ENSMUSG00000001016.13 
            108014596             116282612             116306719             105850014              89427471              90383433 
            
> o <- order(chroms, starts)

> with_introns_grl_test_ordered <- with_introns_grl_test[o]

> head(with_introns_grl_test_ordered)

GRangesList object of length 6:
$ENSMUSG00002075012.1
GRanges object with 2 ranges and 7 metadata columns:
      seqnames          ranges strand |                       source           type     score     phase              gene_id
         <Rle>       <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>          <character>
  [1]     chr3 3069070-3069169      - | dexseq_prepare_annotation.py aggregate_gene        NA      <NA> ENSMUSG00002075012.1
  [2]     chr3 3069070-3069169      - | dexseq_prepare_annotation.py exonic_part           NA      <NA> ENSMUSG00002075012.1
               transcripts exonic_part_number
               <character>        <character>
  [1]                 <NA>               <NA>
  [2] ENSMUST00020182553.1                001
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

$ENSMUSG00002075999.1
GRanges object with 2 ranges and 7 metadata columns:
      seqnames          ranges strand |                       source           type     score     phase              gene_id
         <Rle>       <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>          <character>
  [1]     chr3 3092718-3092817      - | dexseq_prepare_annotation.py aggregate_gene        NA      <NA> ENSMUSG00002075999.1
  [2]     chr3 3092718-3092817      - | dexseq_prepare_annotation.py exonic_part           NA      <NA> ENSMUSG00002075999.1
               transcripts exonic_part_number
               <character>        <character>
  [1]                 <NA>               <NA>
  [2] ENSMUST00020181724.1                001
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

$ENSMUSG00000103416.2
GRanges object with 2 ranges and 7 metadata columns:
      seqnames          ranges strand |                       source           type     score     phase              gene_id
         <Rle>       <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>          <character>
  [1]     chr3 3253093-3253219      + | dexseq_prepare_annotation.py aggregate_gene        NA      <NA> ENSMUSG00000103416.2
  [2]     chr3 3253093-3253219      + | dexseq_prepare_annotation.py exonic_part           NA      <NA> ENSMUSG00000103416.2
               transcripts exonic_part_number
               <character>        <character>
  [1]                 <NA>               <NA>
  [2] ENSMUST00000191986.2                001
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

...
<3 more elements>

> head(with_introns_grl_test)

GRangesList object of length 6:
$ENSMUSG00000000001.5
GRanges object with 18 ranges and 7 metadata columns:
       seqnames              ranges strand |                       source           type     score     phase              gene_id
          <Rle>           <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>          <character>
   [1]     chr3 108014596-108053462      - | dexseq_prepare_annotation.py aggregate_gene        NA      <NA> ENSMUSG00000000001.5
   [2]     chr3 108014596-108016632      - | dexseq_prepare_annotation.py exonic_part           NA      <NA> ENSMUSG00000000001.5
   [3]     chr3 108016633-108016718      - | dexseq_prepare_annotation.py intronic_part         NA      <NA> ENSMUSG00000000001.5
   [4]     chr3 108016719-108016928      - | dexseq_prepare_annotation.py exonic_part           NA      <NA> ENSMUSG00000000001.5
   [5]     chr3 108016929-108019250      - | dexseq_prepare_annotation.py intronic_part         NA      <NA> ENSMUSG00000000001.5
   ...      ...                 ...    ... .                          ...            ...       ...       ...                  ...
  [14]     chr3 108030858-108030999      - | dexseq_prepare_annotation.py  exonic_part          NA      <NA> ENSMUSG00000000001.5
  [15]     chr3 108031000-108031110      - | dexseq_prepare_annotation.py  intronic_part        NA      <NA> ENSMUSG00000000001.5
  [16]     chr3 108031111-108031153      - | dexseq_prepare_annotation.py  exonic_part          NA      <NA> ENSMUSG00000000001.5
  [17]     chr3 108031154-108053203      - | dexseq_prepare_annotation.py  intronic_part        NA      <NA> ENSMUSG00000000001.5
  [18]     chr3 108053204-108053462      - | dexseq_prepare_annotation.py  exonic_part          NA      <NA> ENSMUSG00000000001.5
                transcripts exonic_part_number
                <character>        <character>
   [1]                 <NA>               <NA>
   [2] ENSMUST00000000001.5                001
   [3]                 <NA>                001
   [4] ENSMUST00000000001.5                002
   [5]                 <NA>                002
   ...                  ...                ...
  [14] ENSMUST00000000001.5                007
  [15]                 <NA>                007
  [16] ENSMUST00000000001.5                008
  [17]                 <NA>                008
  [18] ENSMUST00000000001.5                009
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

...
<5 more elements>
```
## Step 4: Generate a GTF file containing exonic and intronic bins.
Now, we would like to convert this GRangesList containing both exonic and intronic bins into a usable GTF file, which can be utilised for downstream read counting. This is done in the final function, `asGFF2`, which takes in a GRangesList and outputs a GTF file. 

### Breakdown of Function 2: asGFF2
```r
> asGFF2 <- function(x) {
  df <- as.data.frame(x) # take in a GRangesList, x
  aggregates <- which(df$type == "aggregate_gene") # returns TRUE for the positions that are aggregate genes
  # character creates a character vector of the specified length
  meta <- character(nrow(df)) # meta is a blank character vector with length = no of rows in df
  # meta has a column called aggregates containing all the gene_id of the aggregates in df
  meta[aggregates] <- sprintf("gene_id \"%s\"", df$gene_id[aggregates])
  # This gives introns a transcript "NA" field, which may not be ideal
  # keeps everything but aggregates
  meta[-aggregates] <- sprintf("transcripts \"%s\"; exonic_part_number \"%s\"; gene_id \"%s\"", df$transcripts[-aggregates], df$exonic_part_number[-aggregates], df$gene_id[-aggregates])
  # paste():Concatenate vectors after converting to character. converts its arguments (via as.character) to character strings
  paste(df$seqnames, "dexseq_prepare_annotation.py", df$type, df$start, df$end, ".", df$strand, ".", meta, sep="\t")
}
```
Breaking down the `asGFF2` function using 1 gene (ENSMUSG00000017688.15):

```r
> df <- as.data.frame(with_introns_grl_test_ordered$ENSMUSG00000017688.15)

> head(df)

  seqnames   start     end  width strand                       source           type score phase               gene_id
1     chr3 3573090 3723112 150023      + dexseq_prepare_annotation.py aggregate_gene    NA    NA ENSMUSG00000017688.15
2     chr3 3573090 3573392    303      + dexseq_prepare_annotation.py    exonic_part    NA    NA ENSMUSG00000017688.15
3     chr3 3573393 3699209 125817      + dexseq_prepare_annotation.py  intronic_part    NA    NA ENSMUSG00000017688.15
4     chr3 3699210 3699407    198      + dexseq_prepare_annotation.py    exonic_part    NA    NA ENSMUSG00000017688.15
5     chr3 3699408 3703118   3711      + dexseq_prepare_annotation.py  intronic_part    NA    NA ENSMUSG00000017688.15
6     chr3 3703119 3703290    172      + dexseq_prepare_annotation.py    exonic_part    NA    NA ENSMUSG00000017688.15
                                transcripts exonic_part_number
1                                      <NA>               <NA>
2                      ENSMUST00000108393.8                001
3                                      <NA>                001
4                      ENSMUST00000108394.3                002
5                                      <NA>                002
6 ENSMUST00000108394.3+ENSMUST00000108393.8                003

> aggregates <- which(df$type == "aggregate_gene") # returns the row numbers of 'aggregate_gene' entries

> head(aggregates)

[1] 1
```

Creating a blank dataframe to subsequently fill with exons:

```r
> meta <- character(nrow(df))

> meta

 [1] "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" ""
```

Fill with aggregate_gene entries first: 

```r
> meta[aggregates] <- sprintf("gene_id \"%s\"", df$gene_id[aggregates])

> meta

 [1] "gene_id \"ENSMUSG00000017688.15\"" ""                                  ""                                 
 [4] ""                                  ""                                  ""                                 
 [7] ""                                  ""                                  ""                                 
[10] ""                                  ""                                  ""                                 
[13] ""                                  ""                                  ""                                 
[16] ""                                  ""                                  ""                                 
[19] ""                                  ""                                  ""                                 
[22] ""                                  ""                                 
```
Populate with the exonic and intronic bins:

```r
> meta[-aggregates] <- sprintf("transcripts \"%s\"; exonic_part_number \"%s\"; gene_id \"%s\"", df$transcripts[-aggregates], df$exonic_part_number[-aggregates], df$gene_id[-aggregates])

> meta

 [1] "gene_id \"ENSMUSG00000017688.15\""                                                                                       
 [2] "transcripts \"ENSMUST00000108393.8\"; exonic_part_number \"001\"; gene_id \"ENSMUSG00000017688.15\""                     
 [3] "transcripts \"NA\"; exonic_part_number \"001\"; gene_id \"ENSMUSG00000017688.15\""                                       
 [4] "transcripts \"ENSMUST00000108394.3\"; exonic_part_number \"002\"; gene_id \"ENSMUSG00000017688.15\""                     
 [5] "transcripts \"NA\"; exonic_part_number \"002\"; gene_id \"ENSMUSG00000017688.15\""                                       
 [6] "transcripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"003\"; gene_id \"ENSMUSG00000017688.15\""
 [7] "transcripts \"NA\"; exonic_part_number \"003\"; gene_id \"ENSMUSG00000017688.15\""                                       
 [8] "transcripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"004\"; gene_id \"ENSMUSG00000017688.15\""
 [9] "transcripts \"NA\"; exonic_part_number \"004\"; gene_id \"ENSMUSG00000017688.15\""                                       
[10] "transcripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"005\"; gene_id \"ENSMUSG00000017688.15\""
[11] "transcripts \"NA\"; exonic_part_number \"005\"; gene_id \"ENSMUSG00000017688.15\""                                       
[12] "transcripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"006\"; gene_id \"ENSMUSG00000017688.15\""
[13] "transcripts \"NA\"; exonic_part_number \"006\"; gene_id \"ENSMUSG00000017688.15\""                                       
[14] "transcripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"007\"; gene_id \"ENSMUSG00000017688.15\""
[15] "transcripts \"NA\"; exonic_part_number \"007\"; gene_id \"ENSMUSG00000017688.15\""                                       
[16] "transcripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"008\"; gene_id \"ENSMUSG00000017688.15\""
[17] "transcripts \"NA\"; exonic_part_number \"008\"; gene_id \"ENSMUSG00000017688.15\""                                       
[18] "transcripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"009\"; gene_id \"ENSMUSG00000017688.15\""
[19] "transcripts \"NA\"; exonic_part_number \"009\"; gene_id \"ENSMUSG00000017688.15\""                                       
[20] "transcripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"010\"; gene_id \"ENSMUSG00000017688.15\""
[21] "transcripts \"NA\"; exonic_part_number \"010\"; gene_id \"ENSMUSG00000017688.15\""                                       
[22] "transcripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"011\"; gene_id \"ENSMUSG00000017688.15\""
[23] "transcripts \"ENSMUST00000108394.3\"; exonic_part_number \"012\"; gene_id \"ENSMUSG00000017688.15\""         
```

Concatenate all the entries together into a larger dataframe containing relevant information:

```r
> paste(df$seqnames, "dexseq_prepare_annotation.py", df$type, df$start, df$end, ".", df$strand, ".", meta, sep="\t")

 [1] "chr3\tdexseq_prepare_annotation.py\taggregate_gene\t3573090\t3723112\t.\t+\t.\tgene_id \"ENSMUSG00000017688.15\""                                                                                    
 [2] "chr3\tdexseq_prepare_annotation.py\texonic_part\t3573090\t3573392\t.\t+\t.\ttranscripts \"ENSMUST00000108393.8\"; exonic_part_number \"001\"; gene_id \"ENSMUSG00000017688.15\""                     
 [3] "chr3\tdexseq_prepare_annotation.py\tintronic_part\t3573393\t3699209\t.\t+\t.\ttranscripts \"NA\"; exonic_part_number \"001\"; gene_id \"ENSMUSG00000017688.15\""                                     
 [4] "chr3\tdexseq_prepare_annotation.py\texonic_part\t3699210\t3699407\t.\t+\t.\ttranscripts \"ENSMUST00000108394.3\"; exonic_part_number \"002\"; gene_id \"ENSMUSG00000017688.15\""                     
 [5] "chr3\tdexseq_prepare_annotation.py\tintronic_part\t3699408\t3703118\t.\t+\t.\ttranscripts \"NA\"; exonic_part_number \"002\"; gene_id \"ENSMUSG00000017688.15\""                                     
 [6] "chr3\tdexseq_prepare_annotation.py\texonic_part\t3703119\t3703290\t.\t+\t.\ttranscripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"003\"; gene_id \"ENSMUSG00000017688.15\""
 [7] "chr3\tdexseq_prepare_annotation.py\tintronic_part\t3703291\t3706282\t.\t+\t.\ttranscripts \"NA\"; exonic_part_number \"003\"; gene_id \"ENSMUSG00000017688.15\""                                     
 [8] "chr3\tdexseq_prepare_annotation.py\texonic_part\t3706283\t3706377\t.\t+\t.\ttranscripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"004\"; gene_id \"ENSMUSG00000017688.15\""
 [9] "chr3\tdexseq_prepare_annotation.py\tintronic_part\t3706378\t3708023\t.\t+\t.\ttranscripts \"NA\"; exonic_part_number \"004\"; gene_id \"ENSMUSG00000017688.15\""                                     
[10] "chr3\tdexseq_prepare_annotation.py\texonic_part\t3708024\t3708130\t.\t+\t.\ttranscripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"005\"; gene_id \"ENSMUSG00000017688.15\""
[11] "chr3\tdexseq_prepare_annotation.py\tintronic_part\t3708131\t3709595\t.\t+\t.\ttranscripts \"NA\"; exonic_part_number \"005\"; gene_id \"ENSMUSG00000017688.15\""                                     
[12] "chr3\tdexseq_prepare_annotation.py\texonic_part\t3709596\t3709751\t.\t+\t.\ttranscripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"006\"; gene_id \"ENSMUSG00000017688.15\""
[13] "chr3\tdexseq_prepare_annotation.py\tintronic_part\t3709752\t3713093\t.\t+\t.\ttranscripts \"NA\"; exonic_part_number \"006\"; gene_id \"ENSMUSG00000017688.15\""                                     
[14] "chr3\tdexseq_prepare_annotation.py\texonic_part\t3713094\t3713181\t.\t+\t.\ttranscripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"007\"; gene_id \"ENSMUSG00000017688.15\""
[15] "chr3\tdexseq_prepare_annotation.py\tintronic_part\t3713182\t3716331\t.\t+\t.\ttranscripts \"NA\"; exonic_part_number \"007\"; gene_id \"ENSMUSG00000017688.15\""                                     
[16] "chr3\tdexseq_prepare_annotation.py\texonic_part\t3716332\t3716484\t.\t+\t.\ttranscripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"008\"; gene_id \"ENSMUSG00000017688.15\""
[17] "chr3\tdexseq_prepare_annotation.py\tintronic_part\t3716485\t3716607\t.\t+\t.\ttranscripts \"NA\"; exonic_part_number \"008\"; gene_id \"ENSMUSG00000017688.15\""                                     
[18] "chr3\tdexseq_prepare_annotation.py\texonic_part\t3716608\t3716844\t.\t+\t.\ttranscripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"009\"; gene_id \"ENSMUSG00000017688.15\""
[19] "chr3\tdexseq_prepare_annotation.py\tintronic_part\t3716845\t3717896\t.\t+\t.\ttranscripts \"NA\"; exonic_part_number \"009\"; gene_id \"ENSMUSG00000017688.15\""                                     
[20] "chr3\tdexseq_prepare_annotation.py\texonic_part\t3717897\t3718019\t.\t+\t.\ttranscripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"010\"; gene_id \"ENSMUSG00000017688.15\""
[21] "chr3\tdexseq_prepare_annotation.py\tintronic_part\t3718020\t3722114\t.\t+\t.\ttranscripts \"NA\"; exonic_part_number \"010\"; gene_id \"ENSMUSG00000017688.15\""                                     
[22] "chr3\tdexseq_prepare_annotation.py\texonic_part\t3722115\t3722904\t.\t+\t.\ttranscripts \"ENSMUST00000108394.3+ENSMUST00000108393.8\"; exonic_part_number \"011\"; gene_id \"ENSMUSG00000017688.15\""
[23] "chr3\tdexseq_prepare_annotation.py\texonic_part\t3722905\t3723112\t.\t+\t.\ttranscripts \"ENSMUST00000108394.3\"; exonic_part_number \"012\"; gene_id \"ENSMUSG00000017688.15\""       
```

We use `lapply` to apply the `asGFF2` function to every gene in the GRangesList.

```r
> final_output <- unlist(lapply(with_introns_grl_test_ordered, asGFF2))

> head(final_output)

                                                                                                                                                   ENSMUSG00002075012.11 
                                                               "chr3\tdexseq_prepare_annotation.py\taggregate_gene\t3069070\t3069169\t.\t-\t.\tgene_id \"ENSMUSG00002075012.1\"" 
                                                                                                                                                   ENSMUSG00002075012.12 
"chr3\tdexseq_prepare_annotation.py\texonic_part\t3069070\t3069169\t.\t-\t.\ttranscripts \"ENSMUST00020182553.1\"; exonic_part_number \"001\"; gene_id \"ENSMUSG00002075012.1\"" 
                                                                                                                                                   ENSMUSG00002075999.11 
                                                               "chr3\tdexseq_prepare_annotation.py\taggregate_gene\t3092718\t3092817\t.\t-\t.\tgene_id \"ENSMUSG00002075999.1\"" 
                                                                                                                                                   ENSMUSG00002075999.12 
"chr3\tdexseq_prepare_annotation.py\texonic_part\t3092718\t3092817\t.\t-\t.\ttranscripts \"ENSMUST00020181724.1\"; exonic_part_number \"001\"; gene_id \"ENSMUSG00002075999.1\"" 
                                                                                                                                                   ENSMUSG00000103416.21 
                                                               "chr3\tdexseq_prepare_annotation.py\taggregate_gene\t3253093\t3253219\t.\t+\t.\tgene_id \"ENSMUSG00000103416.2\"" 
                                                                                                                                                   ENSMUSG00000103416.22 
"chr3\tdexseq_prepare_annotation.py\texonic_part\t3253093\t3253219\t.\t+\t.\ttranscripts \"ENSMUST00000191986.2\"; exonic_part_number \"001\"; gene_id \"ENSMUSG00000103416.2\"" 
```

Lastly, write the entries into a GTF file that can be used. 

```r
> write.table(final_output, file="chr3_flattened.gtf", row.names=F, col.names=F, quote=F)
```

Final GTF file containing intronic bins:

![image](https://user-images.githubusercontent.com/68455070/127977975-f81fb268-9c48-4639-87e0-b94663ee5a32.png)
