#!/usr/bin/env python
import numpy as np
import math
import sys

#lattice constants of super cell is defined by second layer
#limitations
max_strain=1.02
max_atoms=100
max_misangle=0.02
layer_distance=3.25
vacuum_thickness=15
max_multiple=3
min_angle=np.pi/5

if len(sys.argv)==3:
	print_flag=True
elif len(sys.argv)==4:
	print_flag=False
else:
	print('Need 2 or 3 arguments, paths of POSCARs of two unit cells.')
	exit()

def readPOSCAR(path):
	with open(path,'r') as f:
		line=f.read().splitlines()
	car=float(line[1])*np.array([i.split() for i in line[2:5]],dtype='float')
	atom_type=line[5].split()
	atom_num=np.array(line[6].split(),dtype='int')
	if len(set(atom_type))!=len(atom_type): print('Atom types are not correct!');exit()
	atom_dict=dict(zip(atom_type,range(len(atom_type))))
	atom_list=[];start=8
	for val in atom_num:
		end=start+val
		atom_list.append(np.array([i.split() for i in line[start:end]],dtype='float'))
		start=end
	return (car,atom_dict,atom_list)
a_car,a_dict,a_list=readPOSCAR(sys.argv[1])
b_car,b_dict,b_list=readPOSCAR(sys.argv[2])
a=a_car[:2,:2]
b=b_car[:2,:2]
#atoms per cell
numa=sum([len(i) for i in a_list])
numb=sum([len(i) for i in b_list])

cos_angle=np.cos(min_angle)
area_a=abs(np.linalg.det(a)); area_b=abs(np.linalg.det(b))
a1=abs(np.linalg.norm(a[0])); a2=abs(np.linalg.norm(a[1]))
b1=abs(np.linalg.norm(b[0])); b2=abs(np.linalg.norm(b[1]))
cos_a=a[0].dot(a[1])/a1/a2; cos_b=b[0].dot(b[1])/b1/b2
sin_a=np.linalg.det(a)/a1/a2; sin_b=np.linalg.det(b)/b1/b2
b11,b12=b[0];b21,b22=b[1]
max_area=max_atoms*np.square(max_strain)/(numa/area_a+numb/area_b)
max_len=np.sqrt(max_area*max_multiple/np.sin(min_angle))
max_t1=math.ceil(max_len/a1/sin_a)
max_t2=math.ceil(max_len/a2/sin_a)
max_s1=math.ceil(max_len/a1/sin_b)
pset=set()
qset=set()

maxa=np.ceil(max_len); maxb=np.ceil(max_len); maxc=1;
hash_const=10000
def arghash(a,b,c):
	return int((a+maxa)*hash_const)+int((b+maxb)*hash_const)*hash_const*2*maxa+int((c+maxc)*hash_const)*hash_const**2*4*maxa*maxb
result_idx=0;result_list=[]
for t11 in range(0,max_t1+1):
	for t12 in range(max_t2,-max_t2-1,-1):
		p11,p12=t11*a[0]+t12*a[1]; len_p1=np.sqrt(p11*p11+p12*p12)
		if len_p1>max_len or len_p1==0:continue
		for t21 in range(max_t1,-max_t1-1,-1):
			for t22 in range(max_t2,-max_t2-1,-1):
				p21,p22=t21*a[0]+t22*a[1]; len_p2=np.sqrt(p21*p21+p22*p22)
				if (len_p2>max_len) or (len_p2==0) or (len_p2/len_p1>max_multiple) or (len_p1/len_p2>max_multiple) or (p11*p22-p21*p12>max_area) or (p11*p22-p21*p12<0): continue
				cos_p=(p11*p21+p12*p22)/len_p1/len_p2
				if (abs(cos_p)>cos_angle) or (abs(cos_p)>0.01 and cos_p*cos_a<0) or (-0.01<cos_p < 0): continue
				hash_p=arghash(len_p1,len_p2,cos_p)
				if hash_p in pset:
					continue
				else:
					pset.add(hash_p)
				for s11 in range(0,max_s1+1):
					tmp2=2*b11*b12*b21*b22*s11**2 + b22**2*(len_p1**2 - b11**2*s11**2) + b21**2*(len_p1**2 - b12**2*s11**2)
					if tmp2<0: continue
					tmp2=np.sqrt(tmp2)
					tmp1=-b11*b21*s11 - b12*b22*s11
					tmp3=b21*b21+b22*b22
					r1=(tmp1+tmp2)/tmp3; r2=(tmp1-tmp2)/tmp3
					if s11==0: s12_range=set([math.floor(r1),math.ceil(r1)])
					else: s12_range=set([math.floor(r1),math.ceil(r1),math.floor(r2),math.ceil(r2)])
					for s12 in s12_range:
						q11,q12=s11*b[0]+s12*b[1]; len_q1=np.sqrt(q11*q11+q12*q12)
						if len_q1==0 or (len_p1/len_q1>max_strain) or (len_q1/len_p1>max_strain): continue
						cos_pq=(p11*q11+p12*q12)/len_p1/len_q1
						sin_pq=(p11*q12-p12*q11)/len_p1/len_q1
						p_car=np.array([[p11,p12],[p21,p22]])
						R_pq=np.array([[cos_pq,sin_pq],[-sin_pq,cos_pq]])
						pp1,pp2=p_car.dot(R_pq)
						ss21,ss22=np.linalg.inv(b).T.dot(pp2)
						for s21 in set([math.floor(ss21),math.ceil(ss21)]):
							for s22 in set([math.floor(ss22),math.ceil(ss22)]):
								q21,q22=s21*b[0]+s22*b[1]; len_q2=np.sqrt(q21*q21+q22*q22)
								if len_q2==0 or (len_p2/len_q2>max_strain) or (len_q2/len_p2>max_strain): continue
								cos_q=(q11*q21+q12*q22)/len_q1/len_q2
								if abs(cos_p)>1 and abs(cos_p)<1+1e-10: cos_p=np.sign(cos_p)
								if abs(cos_q)>1 and abs(cos_q)<1+1e-10: cos_q=np.sign(cos_q)
								mismatch_angle=abs((np.arccos(cos_p)-np.arccos(cos_q))/np.sin((np.arccos(cos_p)+np.arccos(cos_q))/2))
								if mismatch_angle>max_misangle: continue
								mismatch_len=max(len_p1/len_q1,len_q1/len_p1,len_p2/len_q2,len_q2/len_p2)-1
								atoms_a=int(round(abs((p11*p22-p12*p21)*numa/area_a)))
								atoms_b=int(round(abs((q11*q22-q12*q21)*numb/area_b)))
								q_car=np.array([[q11,q12],[q21,q22]])
								R_q=np.array([[q11,-q12],[q12,q11]])/len_q1
								Q=q_car.dot(R_q)
								#T_p=np.linalg.inv(p_car).dot(q_car).dot(R_q)
								result_list.append([[t11,t12,t21,t22,s11,s12,s21,s22],Q])
								result_idx+=1
								if not print_flag: continue
								print(result_idx,'-------------')
								print("mismatch: len: {:.3f}%, angle: {:.3f}%, total atoms: {}".format(mismatch_len*100,mismatch_angle*100,atoms_a+atoms_b))
								print("{: d} {: d} {: d} {: d}, lattice constant: {:.3f} {:.3f}, angle: {:.3f}, atoms: {}".format(t11,t12,t21,t22,len_p1,len_p2,np.arccos(cos_p)/np.pi*180,atoms_a))
								print("{: d} {: d} {: d} {: d}, lattice constant: {:.3f} {:.3f}, angle: {:.3f}, atoms: {}".format(s11,s12,s21,s22,len_q1,len_q2,np.arccos(cos_q)/np.pi*180,atoms_b))
if len(sys.argv)==4 and int(sys.argv[3])<=len(result_list):
	result_idx=int(sys.argv[3])
else:
	exit()
t11,t12,t21,t22,s11,s12,s21,s22=result_list[result_idx][0]
Q=result_list[result_idx][1]

def zinfo(l):
	sumz=0;numz=0;maxz=-np.inf;minz=np.inf
	for i in l:
		sumz+=i[:,2].sum()
		maxz=max(maxz,i[:,2].max())
		minz=min(minz,i[:,2].min())
		numz+=len(i)
	return (sumz/numz,maxz,minz)
azmean,azmax,azmin=zinfo(a_list)
bzmean,bzmax,bzmin=zinfo(b_list)
lenaz=a_car[2,2];lenbz=b_car[2,2]
lenz=math.ceil(vacuum_thickness+layer_distance+(azmax-azmin)*lenaz+(bzmax-bzmin)*lenbz)
abdist=((azmean-azmin)*lenaz+(bzmax-bzmean)*lenbz+layer_distance)/lenz
Azmean=0.5+abdist/2; Bzmean=0.5-abdist/2
Q=np.array([[Q[0,0],Q[0,1],0],[Q[1,0],Q[1,1],0],[0,0,lenz]])

B_list=[]
S=np.array([[s11,s12,0],[s21,s22,0],[0,0,1]])
for k in range(len(b_list)):
	tmp=[]
	for i in range(min(0,s11,s21,s11+s21),max(0,s11,s21,s11+s21)+1):
		for j in range(min(0,s12,s22,s12+s22),max(0,s12,s22,s12+s22)+1):
			tmp.append(b_list[k]+(i,j,0))
	tmp=np.concatenate(tmp).dot(np.linalg.inv(S))
	tmp=np.round(tmp,decimals=8)
	tmp=tmp[(0<=tmp[:,0]) & (tmp[:,0]<1) & (0<=tmp[:,1]) & (tmp[:,1]<1)]
	tmp[:,2]=(tmp[:,2]-bzmean)*lenbz/lenz+Bzmean
	B_list.append(tmp)
T=np.array([[t11,t12,0],[t21,t22,0],[0,0,1]])
shift_list=np.array([(i,j) for i in [0,0.5] for j in [0,0.5]]).dot(np.linalg.inv(T.dot(S))[:2,:2])
#shift_list=np.array([(i,j) for i in [0] for j in [0]]).dot(np.linalg.inv(T.dot(S))[:2,:2])
for shift_idx,shift in enumerate(shift_list):
	A_list=[]
	for k in range(len(a_list)):
		tmp=[]
		for i in range(min(0,t11,t21,t11+t21),max(0,t11,t21,t11+t21)+1):
			for j in range(min(0,t12,t22,t12+t22),max(0,t12,t22,t12+t22)+1):
				tmp.append(a_list[k]+(i,j,0))
		tmp=np.concatenate(tmp).dot(np.linalg.inv(T))
		tmp=tmp[(0<=tmp[:,0]) & (tmp[:,0]<1) & (0<=tmp[:,1]) & (tmp[:,1]<1)]
		tmp[:,:2]+=shift; tmp[tmp>=1]-=1; tmp[tmp<0]+=1
		tmp[:,2]=(tmp[:,2]-azmean)*lenaz/lenz+Azmean
		A_list.append(tmp)
	A=a_dict.keys();B=b_dict.keys()
	atom_type=[];atom_list=[];atom_num=[]
	for i in A&B:
		atom_type.append(i)
		atom_list.append(A_list[a_dict[i]])
		atom_list.append(B_list[b_dict[i]])
		atom_num.append(len(A_list[a_dict[i]])+len(B_list[b_dict[i]]))
	for i in A-B:
		atom_type.append(i)
		atom_list.append(A_list[a_dict[i]])
		atom_num.append(len(A_list[a_dict[i]]))
	for i in B-A:
		atom_type.append(i)
		atom_list.append(B_list[b_dict[i]])
		atom_num.append(len(B_list[b_dict[i]]))
	atom_list=np.concatenate(atom_list)
	with open('POSCAR_'+str(shift_idx)+'.vasp','w') as f:
		f.write('Hetero\n1\n')
		for i in range(3):
			f.write('   {: 20.16f}  {: 20.16f}  {: 20.16f}\n'.format(*Q[i]))
		for i in atom_type:
			f.write('    {}'.format(i))
		f.write('\n')
		for i in atom_num:
			f.write('    {}'.format(i))
		f.write('\nDirect\n')
		for i in range(len(atom_list)):
			f.write(' {: 18.16f} {: 18.16f} {: 18.16f}\n'.format(*atom_list[i]))
