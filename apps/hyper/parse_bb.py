import sys
import os

def main():
  filename = sys.argv[1]
  name = []
  pre = []
  rank = []
  convert = []
  time = []
  peel = []
  peel_name = []
  with open(filename) as fp:
    append_flag = True
    rank_flag=True
    for line in fp:
      line = line.strip()
      if line == '\n':
        if rank_flag:
          convert.append(-1)
          rank.append(-1)
        append_flag= True
        rank_flag=True
        continue
      split = [x.strip() for x in line.split(':')]
      if split[0] == "preprocess" and append_flag:
        pre.append(float(split[1]))
        append_flag = False
      elif split[0] == "preprocess" and not append_flag:
        pre[len(pre)-1] += float(split[1])
      elif split[0] == "convert":
        convert.append(float(split[1]))
      elif split[0] == "ranking" and rank_flag:
        rank.append(float(split[1]))
        rank_flag = False
      elif split[0] == "ranking" and not rank_flag:
        rank[len(rank)-1] += float(split[1])
      elif split[0] == "Sort" or split[0] == "SortCE" or split[0] == "Hash" or split[0] == "HashCE" or split[0] == "Hist" or split[0] == "HistCE" or split[0] == "Par":
        time.append(float(split[1]))
        name.append(split[0])
      elif split[0] == "E Sort" or split[0] == "E SortCE" or split[0] == "E Hash" or split[0] == "E HashCE" or split[0] == "E Hist" or split[0] == "E HistCE" or split[0] == "E Par":
        time.append(float(split[1]))
        name.append(split[0])
      elif split[0] == "Sort Peel" or split[0] == "Hist Peel" or split[0] == "Par Peel" or split[0] == "Sort E Peel" or split[0] == "Hist E Peel" or split[0] == "Par E Peel":
        peel.append(float(split[1]))
        peel_name.append(split[0])



if __name__ == "__main__":
  main()