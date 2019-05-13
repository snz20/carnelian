# include <stdio.h>
# include <string.h>
# include <stdlib.h>
/*
**  mergeFile
**  Merges existing multiple input files into one output file.
**  
**
**  Return Value:
**  Path of created result file, or 0 in case of bad input or error during merge.
** 
*/

char* mergeFiles(char **fileList, int numFiles, char *outDir, char *fileOutPrefix, char *fileOutExt)
{
    int result = 0;
    FILE *fIn;
    FILE *fOut;
    char buffer[1024 * 1000];
    size_t size;
    size_t read;
    size_t nread;
    size_t written;
    int i,pos;
    char *retVal;

    fOut = NULL;
    result = 1;
    sprintf(buffer, "%s/%s.%s", outDir,fileOutPrefix, fileOutExt);
    fOut = fopen(buffer, "wb");
    size = 0;
    if (fOut == NULL)
    {
    	result *= -1;
    }
    if(result != -1){
        retVal = (char*)malloc(sizeof(char)*(strlen(buffer)+1));
	memcpy(retVal,buffer,strlen(buffer));
	retVal[strlen(buffer)] = '\0';
    	for(i=0;i<numFiles;i++){
	    if(fileList[i] != NULL){
	        fIn = fopen(fileList[i], "rb");
	        if (fIn != NULL)
                {
		    while (!feof(fIn))
                    {
			read = sizeof(buffer);
			read = fread(buffer, 1, read, fIn);
			if (read == 0)
                	{
                    		result *= -1;
                    		break;
                	}
			written = fwrite(buffer, 1, read, fOut);
                	if (written != read)
                	{
                    		result *= -1;
                    		break;
                	}
			size += written;
		    }
	        }
	    }
	    fclose(fIn);
	    if(result == -1)
		break;
	}
	fclose(fOut);
    }
    if(result==-1){
    	return NULL;
    }
    return retVal;
}

int main(int argc, char** argv){
	if(argc<6){
		printf("Insufficient arguments!!!\n./fileMerge outDir outPrefix outExtension inFile1 inFile2 ...\n");
		return -1;
	}
	//mergeFiles(argv[1],argv[2:]);
	//printf("outdir: %s\n",argv[1]);
	//printf("prefix: %s\n",argv[2]);
	//printf("ext: %s\n",argv[3]);
	int count = 0;
	int i,j,k,l;
	while(argv[++count] != NULL);
	//printf("filelist: ");
	char **filelist = (char**)malloc((count-4) * sizeof(char *));
	for(i=4,k=0;i<count;i++,k++){
		l = strlen(argv[i]);
		filelist[k] = (char*)malloc(l+1);
		for(j=0;j<l;j++)
			filelist[k][j] = argv[i][j];
		filelist[k][j]='\0';
		//printf("%s ",filelist[k]);
	}
	//printf("\n");
	char *retval = mergeFiles(filelist, count-4, argv[1], argv[2], argv[3]);
	printf("%s\n",retval);
	return 0;
}

