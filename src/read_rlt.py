#!/usr/bin/python
from scipy import *
class read_data:
   def  __init__(self, data_file,output):
     self.data_file=data_file
     self.output=output
     self.data=[]

   def read_and_ana(self):
     f=open(self.data_file,'r')
     while True:
       line=f.readline()
       if line=='':
         break
       line=line.split()
       if len(line)==0:
         continue
       if line[-1]=="finished" and line[-2]=='0':
         para=line[3:]
         para=self.handle_para(para)
         rlt=f.readline().split()
	 if len(rlt)!=2:
	   continue
         rlt=self.handle_rlt(rlt)
	 self.data.append([para,rlt])
     f.close()
     self.data.sort()
     f=open(self.output,'w')
     for term in self.data:
       #print term
       f.write("%e %e %e %e %e %e %e %e\n"%(term[0][0],term[0][1],\
                term[0][2],term[0][3],term[0][4],term[0][5],term[0][6],term[1]))
     f.close()
   def handle_para(self,para):
     a_l=1.0/double(para[0])
     a_r=1.0/double(para[1])
     c_l=double(para[2])
     c_r=double(para[3])
     b_l=double(para[4])
     b_r=double(para[5])
     r_l=pow(a_l*a_l+b_l*b_l,0.5)
     x_l=pow(r_l,0.5)*pow(2*b_l-r_l,0.5)
     r_r=pow(a_r*a_r+b_r*b_r,0.5)
     x_r=pow(r_r,0.5)*pow(2*b_r-r_r,0.5)
     return [1/a_l,1/a_r,c_l,c_r,b_l,b_r,x_r]
   def handle_rlt(self,rlt):
     return double(rlt[0])
 
for i in xrange(10,22):
  print '%d/rlt.log'%i,'%d/rlt.dat'%i
  a=read_data('%d/rlt.log'%i,'%d/rlt.dat'%i)
  a.read_and_ana()
