#!/usr/bin/python3

#this script will take the TE and gene coordinates to create 50 TE bins of 1 kb distance before and after genes. TEs overlapping with genes are also binned.
#i and last_i are the range iterations for the 50 bins. each iteration has a range of 1kb.
#overlap bins are only written once during the first iteration. Therefore, we include the clause "last_i == 1".
#the "skip" booleans determine if a bin needs to be skipped, if the TEs position doesnt match the bins criteria
#the "loop" booleans determine if TEs are currently summed in a bin. If the new line contains a TE with contradicting criteria to the current bin, the loop is turned off (FALSE), which leads to writing of the last bin into a line of the output. The new line then starts a new bin.
#bin_switch = 1/2/3 determines which bin was written into the output file last. Thereby, the script knows which bins need to be added if bins were skipped. E.g. last bin was pregenic (bin_switch = 1) and the new one will be postgenic: This means we have to add an empty overlap bin (if last_i == 1) or if we have a new gene, we have to add overlap (if last_i == 1) and postgenic of last gene and pregenic and overlap of current gene before writing the line for the postgenic bin of the current gene. Afterwards, the bin_switch will be "3" because last line was postgenic.
#variables with "last_" in front of them are mostly needed for writing of the last bins. The information shouldnt be lossed when a new bin is started. In case there are still missing bins of the last gene, the "last_" variables are necessary to write them down.

with open("gene_coordinates_real_count.txt","r") as f:  #file available as supplementary file in article (Post et al. TBD)
  global lines
  lines = f.readlines()

with open("genome_bin.txt","w") as f:
  f.write(''.join(("caste","\t","no_of_TEs","\t","sum_TE_expression","\t","gene_bin","\t","genome_bin","\t","sum_distance_gene","\t","gene_id","\t","gene_expression","\n")))
  te_exp = []
  distance = []
  gene_id = "Mnat_00001"
  last_gene_id = "Mnat_00001"
  dis_genes = "N/A"
  scaf = "scaffold1"
  last_i = 1
  bin_switch = 3
  pre_skip = post_skip = overlap_skip = pre_loop = post_loop = overlap_loop = False
  last_post_na = last_pre_na = last_overlap_na = False
  pre = post = overlap = na_o = na_pre = na_post = False
  for i in range(1000,50001,1000):
    last_i = i - 999
    print(i)
    print(last_i)
    for line in lines:
      try:
        if (int(line.split()[11])-int(line.split()[10]))/2+int(line.split()[10]) < int(line.split()[8]):
          if last_i == 1:
            for j in range(int(line.split()[10]),int(line.split()[11])+1):
              if j in range(int(line.split()[8]),int(line.split()[9])+1):
                overlap_loop = True
                break
              else:
                post_loop = overlap_loop = False
          else:
            post_loop = overlap_loop = False
        if last_i == 1:
          for j in range(int(line.split()[10]),int(line.split()[11])+1):
            if j in range(int(line.split()[8]),int(line.split()[9])+1):
              pre_loop = post_loop = False
              break
        if (int(line.split()[11])-int(line.split()[10]))/2+int(line.split()[10]) > int(line.split()[9]):
          if last_i == 1:
            for j in range(int(line.split()[10]),int(line.split()[11])+1):
              if j in range(int(line.split()[8]),int(line.split()[9])+1):
                overlap_loop = True
                break
              else:
                pre_loop = overlap_loop = False
          else:
            pre_loop = overlap_loop = False
        if line.split()[3] != gene_id:
          pre_loop = overlap_loop = post_loop = False
        if line.split()[6] or dis_genes != "N/A":
          if dis_genes == "N/A":
            dis_genes = 999999
          if last_gene_id != gene_id:
            last_pre_na = pre_skip
            pre_skip = last_post_na
          if line.split()[7] != scaf:
            post_skip = False
          if pre_loop != True:
            if pre == True:
              pre = False
              if bin_switch == 1:
                if last_i == 1:
                  if last_overlap_na == True:
                    f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                  else:
                    f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                if last_post_na == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              if bin_switch == 2:
                if last_post_na == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              if pre_skip == True:
                f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              else:
                f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              bin_switch = 1
              dis_genes = line.split()[6]
              last_gene_id = gene_id
              last_gene_exp = gene_exp
          if overlap_loop != True:
            if overlap == True:
              overlap = False
              if bin_switch == 1:
                if last_gene_id != gene_id:
                  if last_i == 1:
                    if last_overlap_na == True:
                      f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                    else:
                      f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                  if last_post_na == True:
                    f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                  else:
                    f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                  if pre_skip == True:
                    f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                  else:
                    f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              if bin_switch == 2:
                if last_post_na == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                if pre_skip == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              if bin_switch == 3:
                if pre_skip == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              if last_i == 1:
                if overlap_skip == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              bin_switch = 2
              dis_genes = line.split()[6]
              last_gene_id = gene_id
              last_gene_exp = gene_exp
          if post_loop != True:
            if post == True:
              post = False
              if bin_switch == 1:
                if last_i == 1:
                  if overlap_skip == True:
                    f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                  else:
                    f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              if bin_switch == 2:
                if last_gene_id != gene_id:
                  if last_post_na == True:
                    f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                  else:
                    f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                  if pre_skip == True:
                    f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                  else:
                    f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                  if last_i == 1:
                    if overlap_skip == True:
                      f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                    else:
                      f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              if bin_switch == 3:
                if pre_skip == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                if last_i == 1:
                  if overlap_skip == True:
                    f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                  else:
                    f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              if post_skip == True:
                f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              else:
                f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              dis_genes = line.split()[6]
              last_gene_id = gene_id
              last_gene_exp = gene_exp
              bin_switch = 3
          if line.split()[7] != scaf:
            dis_genes = line.split()[8]
          if line.split()[3] != gene_id:
            last_overlap_na = overlap_skip
            last_post_na = post_skip
            last_pre_na = pre_skip
            overlap_skip = pre_skip = post_skip = False
            try:
              if int(i) > int(line.split()[6]):
                post_skip = True
              else:
                post_skip = False
            except:
              if i > int(dis_genes):
               post_skip = False
          if pre_loop == False and post_loop == False and overlap_loop == False:
            te_exp = []
            distance = []
            pre = post = overlap = False
          caste = line.split()[2]
          gene_id = line.split()[3]
          gene_exp = line.split()[4]
          scaf = line.split()[7]
          gene_start = line.split()[8]
          gene_stop = line.split()[9]
          TE_start = line.split()[10]
          TE_stop = line.split()[11]
          if last_i == 1:
            for j in range(int(line.split()[10]),int(line.split()[11])+1):
              if j in range(int(line.split()[8]),int(line.split()[9])+1):
                caste = line.split()[2]
                gene_id = line.split()[3]
                gene_exp = line.split()[4]
                if float(line.split()[1]) != 0.0:
                  te_exp.append(float(line.split()[1]))
                  distance.append(int(line.split()[5]))
                scaf = line.split()[7]
                gene_start = line.split()[8]
                gene_stop = line.split()[9]
                TE_start = line.split()[10]
                TE_stop = line.split()[11]
                pre_loop = post_loop = False
                overlap_loop = True
                overlap = True
                break
            if overlap == True:
              continue
          if (int(line.split()[11])-int(line.split()[10]))/2+int(line.split()[10]) in [x + 0.5 for x in range(int(line.split()[8])-int(i),int(line.split()[8])-int(last_i)+1)]+[y for y in range(int(line.split()[8])-int(i),int(line.split()[8])-int(last_i)+1)]:
            if pre_skip == True:
              continue
            else:
              caste = line.split()[2]
              gene_id = line.split()[3]
              gene_exp = line.split()[4]
              if float(line.split()[1]) != 0.0:
                te_exp.append(float(line.split()[1]))
                distance.append(int(line.split()[5]))
              scaf = line.split()[7]
              gene_start = line.split()[8]
              gene_stop = line.split()[9]
              TE_start = line.split()[10]
              TE_stop = line.split()[11]
              post_loop = overlap_loop = False
              pre = True
              pre_loop = True
              continue
          if (int(line.split()[11])-int(line.split()[10]))/2+int(line.split()[10]) in [x - 0.5 for x in range(int(line.split()[9])+int(last_i),int(line.split()[9])+int(i)+1)]+[y for y in range(int(line.split()[9])+int(last_i),int(line.split()[9])+int(i)+1)]:
            if post_skip == True:
              continue
            else:
              caste = line.split()[2]
              gene_id = line.split()[3]
              gene_exp = line.split()[4]
              if float(line.split()[1]) != 0.0:
                te_exp.append(float(line.split()[1]))
                distance.append(int(line.split()[5]))
              scaf = line.split()[7]
              gene_start = line.split()[8]
              gene_stop = line.split()[9]
              TE_start = line.split()[10]
              TE_stop = line.split()[11]
              pre_loop = overlap_loop = False
              post = True
              post_loop = True
              continue
          if last_i == 1:
            for j in range(int(line.split()[10]),int(line.split()[11])+1):
              if j in range(int(line.split()[8]),int(line.split()[9])+1):
                overlap = True
                overlap_loop = True
                break
            if overlap == True:
              continue
          if int(line.split()[11]) < int(line.split()[8]) and (int(line.split()[11])-int(line.split()[10]))/2+int(line.split()[10]) not in range(int(line.split()[8]),int(line.split()[9])+1):
            pre = True
            pre_loop = True
            continue
          if int(line.split()[10]) > int(line.split()[9]) and (int(line.split()[11])-int(line.split()[10]))/2+int(line.split()[10]) not in range(int(line.split()[8]),int(line.split()[9])+1):
            post = True
            post_loop = True
            continue
          continue
        if line.split()[3] != gene_id:
          if post == True:
            if bin_switch == 1:
              if last_gene_id != gene_id:
                if last_i == 1:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                if last_post_na == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              if last_pre_na == True:
                f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              else:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              if last_i == 1:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if bin_switch == 2:
              if last_gene_id != gene_id:
                if last_post_na == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                if last_pre_na == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                if last_i == 1:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if bin_switch == 3:
              if last_pre_na == True:
               f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              else:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              if last_i == 1:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if post_skip == True:
              f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            else:
              f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            bin_switch = 3
            last_gene_id = gene_id
            last_gene_exp = gene_exp
            te_exp = []
            last_post_na = last_pre_na =False
            distance = []
          if pre == True:
            if bin_switch == 1:
              if last_gene_id != gene_id:
                if last_i == 1:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                if last_pre_na == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              if last_i == 1:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if bin_switch == 2:
              if last_post_na == True:
                f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              else:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              if last_i == 1:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if last_pre_na == True:
              f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            else:
              f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            bin_switch = 1
            te_exp = []
            last_post_na = last_pre_na = False
            distance = []
            last_gene_id = gene_id
            last_gene_exp = gene_exp
          if overlap == True:
            if bin_switch == 1:
              if last_gene_id != gene_id:
                if last_i == 1:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                if last_post_na == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                if last_pre_na == True:
                  f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                else:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if bin_switch == 2:
              if last_post_na == True:
                f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              else:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              if last_pre_na == True:
                f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              else:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if bin_switch == 3:
              if last_pre_na == True:
                f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              else:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            bin_switch = 2
            te_exp = []
            distance = []
            last_gene_id = gene_id
            last_gene_exp = gene_exp
          if pre == False and post == False and overlap == False:
            if na_o == True or na_pre == True or na_post == True:
              na_o = na_pre = na_post = False
              if last_i == 1:
                if int(line.split()[10]) in range(int(line.split()[8]),int(line.split()[9])+1) and int(line.split()[11]) in range(int(line.split()[8]),int(line.split()[9])+1):
                  caste = line.split()[2]
                  gene_id = line.split()[3]
                  gene_exp = line.split()[4]
                  te_exp.append(float(line.split()[1]))
                  distance.append(int(line.split()[5]))
                  dis_genes = line.split()[6]
                  scaf = line.split()[7]
                  gene_start = line.split()[8]
                  gene_stop = line.split()[9]
                  TE_start = line.split()[10]
                  TE_stop = line.split()[11]
                  overlap = True
                  continue
              if (int(line.split()[11])-int(line.split()[10]))/2+int(line.split()[10]) in [x + 0.5 for x in range(int(line.split()[8])-int(i),int(line.split()[8])-int(last_i)+1)]+[y for y in range(int(line.split()[8])-int(i),int(line.split()[8])-int(last_i)+1)]:
                caste = line.split()[2]
                gene_id = line.split()[3]
                gene_exp = line.split()[4]
                te_exp.append(float(line.split()[1]))
                distance.append(int(line.split()[5]))
                dis_genes = line.split()[6]
                scaf = line.split()[7]
                gene_start = line.split()[8]
                gene_stop = line.split()[9]
                TE_start = line.split()[10]
                TE_stop = line.split()[11]
                pre = True
                continue
              if (int(line.split()[11])-int(line.split()[10]))/2+int(line.split()[10]) in [x - 0.5 for x in range(int(line.split()[9])+int(last_i),int(line.split()[9])+int(i)+1)]+[y for y in range(int(line.split()[9])+int(last_i),int(line.split()[9])+int(i)+1)]:
                caste = line.split()[2]
                gene_id = line.split()[3]
                gene_exp = line.split()[4]
                te_exp.append(float(line.split()[1]))
                distance.append(int(line.split()[5]))
                dis_genes = line.split()[6]
                scaf = line.split()[7]
                gene_start = line.split()[8]
                gene_stop = line.split()[9]
                TE_start = line.split()[10]
                TE_stop = line.split()[11]
                post = True
                continue
          te_exp = []
          distance = []
          pre = post = overlap = False
          caste = line.split()[2]
          gene_id = line.split()[3]
          gene_exp = line.split()[4]
          te_exp.append(float(line.split()[1]))
          distance.append(int(line.split()[5]))
          dis_genes = line.split()[6]
          scaf = line.split()[7]
          gene_start = line.split()[8]
          gene_stop = line.split()[9]
          TE_start = line.split()[10]
          TE_stop = line.split()[11]
        if last_i == 1:
          if int(line.split()[10]) and int(line.split()[11]) in range(int(line.split()[8]),int(line.split()[9])+1):
            if pre == True:
              pre = False
              if bin_switch == 1:
                if last_gene_id != gene_id:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              if bin_switch == 2:
                if last_gene_id != gene_id:
                  f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              bin_switch = 1
              te_exp = []
              last_gene_id = gene_id
              last_gene_exp = gene_exp
              distance = []
              na_o = na_pre = na_post = False
            if post == True:
              post = False
              if bin_switch == 1:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              if bin_switch == 3:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              na_o = na_pre = na_post = False
              bin_switch = 3
              last_gene_id = gene_id
              last_gene_exp = gene_exp
              te_exp = []
              distance = []
            caste = line.split()[2]
            gene_id = line.split()[3]
            gene_exp = line.split()[4]
            te_exp.append(float(line.split()[1]))
            distance.append(int(line.split()[5]))
            dis_genes = line.split()[6]
            scaf = line.split()[7]
            gene_start = line.split()[8]
            gene_stop = line.split()[9]
            TE_start = line.split()[10]
            TE_stop = line.split()[11]
            overlap = True
            continue
        if (int(line.split()[11])-int(line.split()[10]))/2+int(line.split()[10]) in [x + 0.5 for x in range(int(line.split()[8])-int(i),int(line.split()[8])-int(last_i)+1)]+[y for y in range(int(line.split()[8])-int(i),int(line.split()[8])-int(last_i)+1)]:
          if post == True:
            post = False
            if bin_switch == 1:
              if last_gene_id != gene_id:
                if last_i == 1:
                 f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              if last_i == 1:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if bin_switch == 3:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              if last_i == 1:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            na_o = na_pre = na_post = False
            bin_switch = 3
            last_gene_id = gene_id
            last_gene_exp = gene_exp
            te_exp = []
            distance = []
          if overlap == True:
            overlap = False
            if bin_switch == 2:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              if last_na == True:
                f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              else:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if bin_switch == 3:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if last_post_na == True:
              f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              last_post_na = False
              bin_switch = 3
            else:
              bin_switch = 2
            na_o = na_pre = na_post = False
            last_gene_id = gene_id
            last_gene_exp = gene_exp
            te_exp = []
            distance = []
          last_gene_id = gene_id
          last_gene_exp = gene_exp
          caste = line.split()[2]
          gene_id = line.split()[3]
          gene_exp = line.split()[4]
          te_exp.append(float(line.split()[1]))
          distance.append(int(line.split()[5]))
          dis_genes = line.split()[6]
          scaf = line.split()[7]
          gene_start = line.split()[8]
          gene_stop = line.split()[9]
          TE_start = line.split()[10]
          TE_stop = line.split()[11]
          pre = True
          continue
        if (int(line.split()[11])-int(line.split()[10]))/2+int(line.split()[10]) in [x + 0.5 for x in range(int(line.split()[9])+int(last_i),int(line.split()[9])+int(i)+1)]+[y for y in range(int(line.split()[9])+int(last_i),int(line.split()[9])+int(i)+1)]:
          if pre == True:
            pre = False
            if bin_switch == 1:
              if last_i == 1:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
            if bin_switch == 2:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
            f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            na_o = na_pre = na_post = False
            bin_switch = 1
            last_gene_id = gene_id
            last_gene_exp = gene_exp
            te_exp = []
            distance = []
          if overlap == True:
            overlap = False
            if bin_switch == 2:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              if last_na == True:
                f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              else:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if bin_switch == 3:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if last_post_na == True:
              f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              last_post_na = False
              bin_switch = 3
            else:
              bin_switch = 2
            last_gene_id = gene_id
            last_gene_exp = gene_exp
            na_o = na_pre = na_post = False
            te_exp = []
            distance = []
          last_gene_id = gene_id
          last_gene_exp = gene_exp
          caste = line.split()[2]
          gene_id = line.split()[3]
          gene_exp = line.split()[4]
          te_exp.append(float(line.split()[1]))
          distance.append(int(line.split()[5]))
          dis_genes = line.split()[6]
          scaf = line.split()[7]
          gene_start = line.split()[8]
          gene_stop = line.split()[9]
          TE_start = line.split()[10]
          TE_stop = line.split()[11]
          post = True
          continue
        else:
          continue
      except:
        if pre == True:
          if bin_switch == 1:
            if last_i == 1:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
            if last_post_na == True:
              f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
            else:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
          if bin_switch == 2:
            if last_post_na == True:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
            else:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
          if last_post_na == True:
            f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
          else:
            f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
          last_post_na = last_pre_na = na_o = na_pre = na_post = False
          bin_switch = 1
          last_gene_id = gene_id
          last_gene_exp = gene_exp
        if overlap == True:
          if bin_switch == 2:
            if last_post_na == True:
              f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
            else:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
            if last_pre_na == True:
              f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            else:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
          if bin_switch == 3:
            if last_pre_na == True:
              f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            else:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
          f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
          bin_switch = 2
          last_post_na = last_pre_na = na_o = na_pre = na_post = False
          last_gene_id = gene_id
          last_gene_exp = gene_exp
        if post == True:
          if bin_switch == 1:
            if last_gene_id != gene_id:
              if last_i == 1:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              if last_post_na == True:
                f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              else:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(last_gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("0"),"\t",str(last_gene_id),"\t",str(last_gene_exp),"\n")))
              if last_pre_na == True:
                f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
              else:
                f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if last_i == 1:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
          if bin_switch == 3:
            if last_pre_na == True:
              f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            else:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_pre_"+str(i),"\t",str(str("-")+str(last_i)+str(" to ")+str("-")+str(i)),"\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
            if last_i == 1:
              f.write(''.join((str(caste),"\t",str("0"),"\t",str("0"),"\t",str(gene_id)+"_overlap","\t","overlap","\t",str("0"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
          if last_post_na == True:
            f.write(''.join((str(caste),"\t",str("N/A"),"\t",str("N/A"),"\t",str(gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str("N/A"),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
          else:
            f.write(''.join((str(caste),"\t",str(len(te_exp)),"\t",str(sum(te_exp)),"\t",str(gene_id)+"_post_"+str(i),"\t",str(str(last_i)+str(" to ")+str(i)),"\t",str(sum(distance)),"\t",str(gene_id),"\t",str(gene_exp),"\n")))
          last_post_na = last_pre_na = na_o = na_pre = na_post = False
          bin_switch = 3
          last_gene_id = gene_id
          last_gene_exp = gene_exp
        te_exp = []
        distance = []
        pre = post = overlap = False
