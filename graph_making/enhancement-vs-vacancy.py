import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


#Times New Roman
tnr= fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman.ttf')
#Times New Roman italic
tnri = fm.FontProperties(fname='/home/pxiang/extrafonts/Times_New_Roman_Italic.ttf')

average=np.array([])
percentage, enhancement, enhancement2 = np.loadtxt("enhancement_target_time.txt",usecols=(0,4,5)).T
num_entry,  = np.shape(enhancement) 
print num_entry

mapping = {2.5:'a1', 5:'a2', 7.5:'a3', 10:'a4', \
           12.5:'a5', 15:'a6', 17.5:'a7', 20:'a8', \
           30:'a9', 40:'a10', 50:'a11', 60:'a12', \
           70:'a13', 80:'a14', 90:'a15'}
# 2.5
a1=[] #np.array([])
# 5
a2=[] #np.array([])
# 7.5
a3=[] #np.array([])
# 10
a4=[] #np.array([])
# 12.5
a5=[] #np.array([])
# 15
a6=[] #np.array([])
# 17.5
a7=[] #np.array([])
# 20
a8=[] #np.array([])

a9=[]
a10=[]
a11=[]
a12=[]
a13=[]
a14=[]
a15=[]


  
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

a1=[] #np.array([])
# 5
a2=[] #np.array([])
# 7.5
a3=[] #np.array([])
# 10
a4=[] #np.array([])
# 12.5
a5=[] #np.array([])
# 15
a6=[] #np.array([])
# 17.5
a7=[] #np.array([])
# 20
a8=[] #np.array([])

a9=[]
a10=[]
a11=[]
a12=[]
a13=[]
a14=[]
a15=[]

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
plt.errorbar(percentage_reduced, average_enhancement, yerr=errors, fmt='bo', capsize=4, markersize=3)
# to shift the plot a little bit
percentage_reduced2 = np.array(percentage_reduced2)
percentage_reduced2 = percentage_reduced2 #+ 0.5
plt.errorbar(percentage_reduced2, average_enhancement2, yerr=errors2, fmt='ro', capsize=4, markersize=3)
#plt.xlim(-0.7, 95)
plt.xlim(-10,95)
plt.ylim(-10, 80)
plt.xticks([0,10,20,30,40,50,60,70,80,90], fontsize=20, fontproperties=tnr)
plt.yticks([1,10,20,30,40,50,60,70,80], fontsize=20, fontproperties=tnr)
# plot a dot
# t_target = 290, p_290.dat

command='''cd /home/pxiang/run_jobs/wavepacket_1d/wavepacket_focusing_2d/Ja_22.83_kHz/time_0.5ms_sigma_1d-3_no_focusing;
            sed -n '/ 0 *0 /p' p_290.dat | awk '{ print $3 }' '''
p_no_focusing=os.popen(command).read()
p_no_focusing=float(p_no_focusing)

command='''cd /home/pxiang/run_jobs/wavepacket_1d/wavepacket_focusing_2d/Ja_22.83_kHz/time_0.5ms_sigma_1d-3_focus_51.51;
            sed -n '/ 0 *0 /p' p_290.dat | awk '{ print $3 }' '''
p_focusing=os.popen(command).read()
p_focusing=float(p_focusing)

command='''cd /home/pxiang/run_jobs/wavepacket_1d/wavepacket_focusing_2d/Ja_22.83_kHz/time_0.5ms_sigma_1d-3_focus_51.51;
            sed -n '/ 0 *0 /p' p_0.dat | awk '{ print $3 }' '''
p_initial=os.popen(command).read() 
p_initial=float(p_initial)

# P_focusing(t_target)/P_no_focusing(t_target)                     
plt.plot(0, p_focusing/p_no_focusing, 'bo')
# P_focusing(t_target)/P(t=0)
plt.plot(0, p_focusing/p_initial, 'ro')
# plot a line with heigh 1 
plt.axhline(y=1,linestyle='-',color='black', linewidth=0.05)
#plt.text(-5, 2.0, "1", color='green', fontsize=18, fontproperties=tnr)
plt.xlabel('Vacancy percentage (%)', fontsize=24, fontproperties=tnr)
plt.ylabel('Enhancement factor', fontsize=24, fontproperties=tnr)

# anotation
#plt.text(15,75, r"$T_{target}$ is the time when the probability $P$"+"\n at the target molecule (51,51) is maximal" + "\n when there is no vacancy", ha='left', va='top')
#plt.text(15,60,r'blue $\rightarrow$ $P_{focusing}(T_{target})/P_{no \;focusing}(T_{target})$', ha='left', va='top',color='blue',fontsize=19) 
#plt.text(15,50,r'red $\rightarrow$ $P_{focusing}(T_{target})/P(t=0)$', ha='left', va='top',color='red',fontsize=19) 
#plt.text(37,37, 'Note the red dots have been shift a little \n to the right for better visualization')       
plt.savefig("enhancement-vs-vacancy.pdf")
plt.show()
