#special entries for the dictionnary passed by :
#'directory' is the directory where to put file, if not found, it is set to '../Data/'
#'startName' is the file startName, if not found, this is set to  'SingleMu_upscope'
#'extension' is the file extension, if not found, this is set to '.root'
def generateFileName(paramsBydict):
  filename='../Data/'
  if ('directory' in paramsBydict): 
    filename=paramsBydict['directory']
  if ('startName' in paramsBydict):
    filename=filename+paramsBydict['startName']
  else:
    filename=filename+'SingleMu_upscope'
  sortedKeys=sorted(paramsBydict)
  for x in sortedKeys:
    if (x not in ['directory','startName','nevents','extension']):
       filename=filename+'_'+x+str(paramsBydict[x])
  if ('nevents' in paramsBydict):
    filename=filename+'_'+str(paramsBydict['nevents'])
  if ('extension' in paramsBydict):
    filename=filename+str(paramsBydict['extension'])
  else:
    filename=filename+'.root'
  return filename


if __name__ == "__main__":
   a=dict()
   a['Pt']=60
   a['zvtx']=30
   a['etamin']=2.3
   a['nevents']=1000
   print generateFileName(a)

