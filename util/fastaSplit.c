# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <sys/stat.h>
# include <sys/types.h>
# include <unistd.h>
/*
**  splitFile
**  Splits an existing input file into multiple output files with a specified
**  maximum file size.
**
**  Return Value:
**  Number of created result files, or 0 in case of bad input data or a negative
**  value in case of an error during file splitting.
*/
char *strrev(char *str)
{
      char *p1, *p2;

      if (! str || ! *str)
            return str;
      for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
      {
            *p1 ^= *p2;
            *p2 ^= *p1;
            *p1 ^= *p2;
      }
      return str;
}

int safeCreateDir(char *parent, int num){
	char dirname[10000];
	sprintf(dirname,"%s/split%03d", parent, num);
	struct stat st = {0};
	if (stat(dirname, &st) == -1){
		mkdir(dirname,0755);
		//printf("%s: ", dirname);
		return 0;
	}
	//printf("%s directory already exists!!\n", dirname);
	return 0;
}

int splitFile(char *fileIn, char *outDir, char *prefixOut, size_t maxSize)
{
    int result = 0;
    FILE *fIn;
    FILE *fOut;
    char buffer[1024 * 1000], rem[10000];
    size_t size;
    size_t read;
    size_t nread;
    size_t written;
    int pos, dirres;
    rem[0] = '\0';
    if ((fileIn != NULL) && (maxSize > 0))
    {
        fIn = fopen(fileIn, "rb");
        if (fIn != NULL)
        {
            fOut = NULL;
            //result = 1;   /* we have at least one part */

            while (!feof(fIn))
            {
                /* initialize (next) output file if no output file opened */
                if (fOut == NULL)
                {
		    dirres = safeCreateDir(outDir, result);
		    if(dirres != -1){
                    	sprintf(buffer, "%s/split%03d/%s.%03d.fasta", outDir, result, prefixOut, result);
                    	fOut = fopen(buffer, "wb");
			//printf("Created %d\n",result);
			result++;
			if(strlen(rem)>0){
                        	written = fwrite(rem, 1, strlen(rem), fOut);
                        	size = strlen(rem);
				rem[0] = '\0';
                    	}
                    }
		    if (fOut == NULL)
                    {
                        result *= -1;
                        break;
                    }
                    else{
			size = 0;
		    }
                }
		//printf("size: %zu\n",size);
                /* calculate size of data to be read from input file in order to not exceed maxSize */
                read = sizeof(buffer);
                if ((size + read) > maxSize)
                {
                    read = maxSize - size;
                }
		//printf("read: %zu\n",read);
                /* read data from input file */
                nread = fread(buffer, 1, read, fIn);
		//printf("nread: %zu\n",nread);
                if (read == 0)
                {
                    result *= -1;
                    break;
                }
		//printf("read: %zu\nnread: %zu\nbefore: %c\n",read,nread,buffer[nread-1]);
		/* additional check for fasta record end */
		pos=0;
		buffer[nread] = '\0';
		//printf("buffer:%s\n",buffer);
		//printf("buffer length: %zu\n",strlen(buffer));
		//rem = (char *) malloc(sizeof(char));
		if(nread == read){
			while(buffer[nread-1] != '>'){
				//realloc(rem, (sizeof(char)));
				rem[pos++] = buffer[nread-1];
				nread--;
			}
			//printf("middle: read: %zu\nnread: %zu\nafter: %c\n",read,nread,buffer[nread-1]);
			if(buffer[nread-1] == '>'){
				nread--;
			}
			rem[pos++]='>';
			rem[pos]='\0';
			//printf("%s\n",rem);
			strrev(rem);
			//printf("pass\n*****\n%s\n",rem);
			//printf("read: %zu\nnread: %zu\nafter: %c\n",read,nread,buffer[nread-1]);
                	/* write data to output file */
		}
                written = fwrite(buffer, 1, nread, fOut);
                if (written != nread)
                {
                    result *= -1;
                    break;
                }

                /* update size counter of current output file */
                size += written;
		//printf("size: %zu\n",size);
                if (size >= maxSize-(read-nread))   /* next split? */
                {
		    //printf("has %zu bytes\n",size);
                    fclose(fOut);
                    fOut = NULL;
                    //result++;
		    //sprintf(buffer, "%s.%03d", fileIn, result);
		    //printf("result: %d\n", result);
		    //dirres = safeCreateDir(outDir, result-1);
                    //if(dirres != -1){
                    //    sprintf(buffer, "%s/split%03d/%s.%03d.fasta", outDir, result-1, prefixOut, result-1);
                    //    fOut = fopen(buffer, "wb");
                    //}
                    //if (fOut == NULL)
                    //{
                    //    result *= -1;
                    //    break;
                    //}
		    //printf("being written:\n********\n%s\n",rem);
		    //printf("%s ",buffer);
		    //written = fwrite(rem, 1, strlen(rem), fOut);
		    //size = strlen(rem);
                }
		else{
			written = fwrite(rem, 1, strlen(rem), fOut);
			size += written;
			//printf("being written:\n********\n%s\n",rem);
			//printf("rem length: %zu\n", strlen(rem));
			rem[0] = '\0';
			//printf("rem length after: %zu\n", strlen(rem));
		}
            }

            /* clean up */
            if (fOut != NULL)
            {
		//printf("has %zu bytes\n",size);
                fclose(fOut);
            }
            fclose(fIn);
        }
    }

    return (result);
}

int main(int argc, char** argv){
	if(argc<5){
		printf("Insufficient arguments!!\nUsage: ./fastaSplit inFasta outDir outPrefix maxSize\n");
		return -1;
	}
	int maxSize = atoi(argv[4]);
	int result = splitFile(argv[1],argv[2],argv[3],maxSize);
	printf("%d\n",result);
	return 0;
}

