# Copyright (c) 2024 Mitsuru Ohno
# Use of this source code is governed by a BSD-3-style
# (license that can be found in the LICENSE file. )

# 08/13/2024, M. Ohno
#Markush構造から、SMARTS を発生させる構造発生器の機能を有する、Pythonの関数
#NOTE:  
#1. 改変したい置換基を"Q<一桁添数>"で記載した骨格のSMARTSを記述する  
#2. 全ての可変置換基で、候補官能基が同一な場合には、それをリストとして与える。異なる場合には、可変置換基をキー、その候補官能基のリストを値とする辞書で与える  
#可変の時間機をQ<1桁の添数>で定義する。Qのみや、Q<2桁以上の添数>は不可で、エラーとなる。よって、可変の置換機は添数0-9の10種類まで
#置換機に水素を含む場合には[H]で指定する。

# A Python function as a structure generator to generate SMARTS/SMILES from Markush structures.
#NOTE:
#Define a variable timer with Q<1-digit index>.
#Only Q or Q<addition of 2 or more digits> is not allowed and will result in an error.
#Therefore, there are up to 10 types of variable replacement machines with indexes 0-9.
#If the displacement machine contains hydrogen, specify with [H].

import sys
import re
import itertools

def markushy(skelton, Q):
  Qn = re.findall('Q[0-9]{0,100}', skelton)
  for e in Qn: #可変置換基の添数チェック
    if len(e) != 2: #添数を2桁化するのであれば、この長さを調整
      print('Please review the definition of Q (index).')
      #return
      sys.exit()
  if type(Q) == list: #全ての可変置換基の候補置換基が同一の場合
    Qdic = {k:Q for k in Qn}
  elif type(Q) == dict: #可変置換基の候補官能基が異なる場合
    if set(Q.keys()) != set(Qn):
      print('Please review the definition of Q (dictionary key).')
      #return
      sys.exit()
    else:
      Qdic = Q
  else:
    print('Please review the definition of Q (dictionary).')
    #return
    sys.exit()
  Qcomb = [fg for fg in itertools.product(*Qdic.values())] #可変置換基同志の候補官能基の全ての組み合わせを発生
  Qcomb_dic = [{k:v for k, v in zip(Qn, e)} for e in Qcomb]
  rev_skelton = []
  for d in Qcomb_dic:
    inter_skl = re.sub("({})".format("|".join(map(re.escape, Qcomb_dic[0].keys()))), lambda inter_skl: d[inter_skl.group()], skelton)
    #SMARTS/SMILES中の可変置換基を単語として検索し、候補官能基文字列（単語）と置換

  #[H]をsome atome with some number of hydrogen (ex. [CX3H2]) 形式に修正する
    while re.search('\(?\[H\]\)?', inter_skl):
      #文字列を逆転して、[H]の前の原子を探す
      inv_inter_skl = inter_skl[::-1]
      ref_position = len(inter_skl)-1 #lenの戻り値は1始まりなので0始まりに修正
      ref_point = re.search('\(?\[H\]\)?', inter_skl) #[H]記述の始点終点
      inv_start = ref_position - ref_point.end()+1 #逆転文字列での[H]記述の始点終点
      inv_end = ref_position - ref_point.start()
      inv_range_start = inv_inter_skl.find(']', inv_end) #[H]記述直前の原子の記述範囲
      inv_range_end = inv_inter_skl.find('[', inv_end+1)


      targ_range_start = ref_position-inv_range_end #変更する原子の記述の始点終点
      targ_range_end = ref_position-inv_range_start
      form_targ = inter_skl[:targ_range_start] #修正する原子より前のSMARTS文字列
      follow_targ = inter_skl[ref_point.end():] #修正[H]より後ろのSMARTS文字列
      targ_strings = inter_skl[targ_range_start:targ_range_end+1] #修正対象の原子

      #Hの数の指定がある場合、各1を増やす
      if re.search('H[0-9]' , targ_strings):
        targ_strings = re.sub('H3', 'H4', targ_strings) #H多い方から修正
        targ_strings = re.sub('H2', 'H3', targ_strings)
        targ_strings = re.sub('H1', 'H2', targ_strings)
        targ_strings = re.sub('H0', 'H1', targ_strings)
      #絞り込んだ原子にHの数指定がない場合、H1を追記する
      else:
        atoms = '[B|C|N|O|Al|Si|P|S|Ge|As|Se]X?[0-9]?' #共有結合なので典型元素非金属
        atm_patt = re.compile(atoms)
        atm_l = atm_patt.findall(targ_strings)
        for e in atm_l:
          targ_strings = re.sub(e, e+'H1', targ_strings)

      inter_skl = form_targ + targ_strings + follow_targ

    rev_skelton.append(inter_skl)
  return rev_skelton
