import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


#Times New Roman
tnr= fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman.ttf')
#Times New Roman italic
tnri = fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman_Italic.ttf')

average=np.array([])
percentage, enhancement2, enhancement = np.loadtxt("enhancement_t_matrix.txt",usecols=(0,2,3)).T
num_entry,  = np.shape(enhancement) 
print num_entry

mapping = {10:'a1', 20:'a2', 30:'a3', 40:'a4', \
           50:'a5', 60:'a6', 70:'a7', 80:'a8', \
           90:'a9'}
# 10 percent
a1=[] #np.array([])

# 20 percent
a2=[] #np.array([])

# 30 percent
a3=[] #np.array([])

# 40 percent
a4=[] #np.array([])

# 50 percent
a5=[] #np.array([])

# 60 percent
a6=[] #np.array([])

# 70 percent
a7=[] #np.array([])

# 80 percent
a8=[] #np.array([])

# 90 percent
a9=[]



  
for i in np.arange(num_entry):
  array_name = mapping[percentage[i]]
  tmp=eval(array_name)
  tmp.append(enhancement[i]) 
  

percentage_reduced=[]
average_enhancement=[]
errors = []

for perc in mapping.keys():
  tmp = eval(mapping[perc])
  tmp = np.asarray(tmp)
  # get rid of the largest and smallest values
  index=[np.argmax(tmp), np.argmin(tmp)] 
  tmporary=np.delete(tmp, index)
  tmp=tmporary
  
  n, = np.shape(tmp)
  print perc, n
  percentage_reduced.append(perc)
  average_enhancement.append(np.sum(tmp)/n)
  sd = np.std(tmp,ddof=1) # standard deviation
  errors.append(1.96*sd/np.sqrt(n*1.0))
#print type(a1)

# 10 percent
a1=[] #np.array([])

# 20 percent
a2=[] #np.array([])

# 30 percent
a3=[] #np.array([])

# 40 percent
a4=[] #np.array([])

# 50 percent
a5=[] #np.array([])

# 60 percent
a6=[] #np.array([])

# 70 percent
a7=[] #np.array([])

# 80 percent
a8=[] #np.array([])

# 90 percent
a9=[]


for i in np.arange(num_entry):
  array_name = mapping[percentage[i]]
  tmp2=eval(array_name)
  tmp2.append(enhancement2[i]) 

percentage_reduced2=[]
average_enhancement2=[]
errors2 = []

for perc2 in mapping.keys():
  tmp2 = eval(mapping[perc2])
  tmp2 = np.asarray(tmp2)
  # get rid of the largest and smallest values
  index=[np.argmax(tmp2), np.argmin(tmp2)] 
  tmporary=np.delete(tmp2, index)
  tmp2=tmporary 
  
  n2, = np.shape(tmp2)
  print perc2, n2
  percentage_reduced2.append(perc2)
  average_enhancement2.append(np.sum(tmp2)/n2)
  sd2 = np.std(tmp2,ddof=1) # standard deviation
  errors2.append(1.96*sd2/np.sqrt(n2*1.0))

np.savetxt("data.txt", np.transpose((percentage_reduced, average_enhancement, errors)))
plt.figure()
plt.subplots_adjust(bottom=0.13,top=0.93)
plt.tick_params(direction='out', top='off', right='off',length=8,width=1.2,pad=8)
#plt.errorbar(percentage_reduced, average_enhancement, yerr=errors, fmt='bo', capsize=4, markersize=3)
# to shift the plot a little bit
percentage_reduced2 = np.array(percentage_reduced2)
percentage_reduced2 = percentage_reduced2 #+ 0.5
plt.errorbar(percentage_reduced2, average_enhancement2, yerr=errors2, fmt='ro', capsize=4, markersize=3)

plt.xlim(-10,100)
plt.ylim(-10, 350)
plt.xticks([0,10,20,30,40,50,60,70,80,90,100], fontsize=20, fontproperties=tnr)
plt.yticks([1,50,100,150,200,250,300,350], fontsize=20, fontproperties=tnr)
plt.yticks(fontsize=20, fontproperties=tnr)
# plot a dot
# t_target = 290, p_290.dat

#command='''cd /home/pxiang/run_jobs/wavepacket_1d/wavepacket_focusing_2d/T_matrix_focusing/with_focusing/0_percent;
            #sed -n '/ 0 *0 /p' p_4000_no.dat | awk '{ print $3 }' '''
#p_no_focusing=os.popen(command).read()
#p_no_focusing=float(p_no_focusing)

command='''cd /home/pxiang/run_jobs/wavepacket_1d/wavepacket_focusing_2d/T_matrix_focusing/with_focusing/0_percent;
            sed -n '/ 0 *0 /p' p_4000.dat | awk '{ print $3 }' '''
p_focusing=os.popen(command).read()
p_focusing=float(p_focusing)

command='''cd /home/pxiang/run_jobs/wavepacket_1d/wavepacket_focusing_2d/T_matrix_focusing/with_focusing/0_percent;
            sed -n '/ 0 *0 /p' p_0.dat | awk '{ print $3 }' '''
p_initial=os.popen(command).read() 
p_initial=float(p_initial)

# P_focusing(t_target)/P_no_focusing(t_target)                     
#plt.plot(0, p_focusing/p_no_focusing, 'bo')
# P_focusing(t_target)/P(t=0)
plt.plot(0, p_focusing/p_initial, 'ro')
## plot a line with heigh 1 
plt.axhline(y=1,linestyle='-',color='black', linewidth=0.05)
#plt.text(-5, 2.0, "1", color='green', fontsize=18, fontproperties=tnr)
plt.xlabel('Vacancy percentage (%)', fontsize=24, fontproperties=tnr)
plt.ylabel('Enhancement factor', fontsize=24, fontproperties=tnr)

# anotation
#plt.text(15,75, r"$T_{target}$ is the time when the probability $P$"+"\n at the target molecule (51,51) is maximal" + "\n when there is no vacancy", ha='left', va='top')
#plt.text(15,60,r'blue $\rightarrow$ $P_{focusing}(T_{target})/P_{no \;focusing}(T_{target})$', ha='left', va='top',color='blue',fontsize=19) 
#plt.text(15,50,r'red $\rightarrow$ $P_{focusing}(T_{target})/P(t=0)$', ha='left', va='top',color='red',fontsize=19) 
#plt.text(37,37, 'Note the red dots have been shift a little \n to the right for better visualization') 
plt.savefig("enhancement-vs-vacancy-t-matrix-other.png")      
plt.savefig("enhancement-vs-vacancy-t-matrix-other.pdf")
plt.show()
