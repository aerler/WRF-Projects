import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.font_manager as font_manager


t31nc = Dataset("./PHIS_T31.nc", 'r')
t85nc = Dataset("./PHIS_T85.nc", 'r')
v11nc = Dataset("./PHIS_0.9x1.25.nc", 'r')
etoponc = Dataset("./topo10min.nc", 'r')

zt31 = t31nc.variables["PHIS"][0,:,:]/9.81
latt31 = t31nc.variables["lat"][:]
lont31 = t31nc.variables["lon"][:]

zt85 = t85nc.variables["PHIS"][0,:,:]/9.81
latt85 = t85nc.variables["lat"][:]
lont85 = t85nc.variables["lon"][:]

zv11 = v11nc.variables["PHIS"][0,:,:]/9.81
latv11 = v11nc.variables["lat"][:]
lonv11 = v11nc.variables["lon"][:]

zetopo = etoponc.variables["htopo"][:,:]
latetopo = etoponc.variables["lat"][:]
lonetopo = etoponc.variables["lon"][:]

fig = plt.figure(facecolor='white', figsize = (11, 8.5)) #landscape
#fig = plt.figure(facecolor='white', figsize = (8.5,11)) #portrait

ax = fig.add_subplot(111)

print zt31.shape
print zt85.shape
print zv11.shape


#labels = ['T31 at 75 N','T85 at 75 N','FV1x1 at 75 N','T31 at 68 N','T85 at 68 N','FV1x1 at 68 N']
labels = ['T31 at 75 N','T85 at 75 N','FV1x1 at 75 N','ETOPO 10-minute at 75 N','T31 at 68 N','T85 at 68 N','FV1x1 at 68 N']


lat75t31 = np.where(latt31[:] >= 75.0)[0][0]
lat75t85 = np.where(latt85[:] >= 75.0)[0][0]
lat75v11 = np.where(latv11[:] >= 75.0)[0][0]
lat75etopo = np.where(latetopo[:] >= 75.0)[0][0]
lat68t31 = np.where(latt31[:] >= 68.0)[0][0]
lat68t85 = np.where(latt85[:] >= 68.0)[0][0]
lat68v11 = np.where(latv11[:] >= 68.0)[0][0]
lat68etopo = np.where(latetopo[:] >= 68.0)[0][0]

alp=1.0
lines = ax.plot(lont31,zt31[lat75t31,:],linewidth=2,color='red',label=labels[0],alpha=alp)
lines.append(ax.plot(lont85,zt85[lat75t85,:],linewidth=2,color='blue',label=labels[1],alpha=alp))
lines.append(ax.plot(lonv11,zv11[lat75v11,:],linewidth=2,color='green',label=labels[2],alpha=alp))
lines.append(ax.plot(lonetopo,zetopo[lat75etopo,:],linewidth=2,color='black',label=labels[3],alpha=alp))

#lines.append(ax.plot(lont31,zt31[lat68t31,:],linewidth=2,color='orange',label=labels[3],alpha=alp))
#lines.append(ax.plot(lont85,zt85[lat68t85,:],linewidth=2,color='lightblue',label=labels[4],alpha=alp))
#lines.append(ax.plot(lonv11,zv11[lat68v11,:],linewidth=2,color='lightgreen',label=labels[5],alpha=alp))


ax.set_xlim([180.0,360.0])

ax.set_ylabel('Height (m)', color='k')
#ax.yaxis.set_label_coords(-0.09, 0.23)
#axmoc.yaxis.set_label_coords(-0.115, 0.77)

titlestr="Atmospheric Model Topographic Height"
plt.title(titlestr)

ax.set_xlabel('Longitude')

props = font_manager.FontProperties(size=20)
leg2 = ax.legend(loc='upper left', shadow=True, fancybox=True, prop=props)
leg2.get_frame().set_alpha(0.5)


plt.show()
