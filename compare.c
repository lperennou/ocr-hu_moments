#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <png.h>
#include <errno.h>
#include <getopt.h>



static struct option long_options[]={
   {"help",no_argument,0,'h'},
   {"version",no_argument,0,'v'},
   {0,0,0,0}
};

#define VERSION "Final version written by Mathilde, Kevin and Loïc with Mr. Pujol's help\n"

void Version() {
   printf("\n%s\n\n",VERSION);
}
void Help() {
   printf("\nProvide a file name or a file name and database name as argument. The program will give you the 6 first Hu invariant moments and the closest figure from your database\n \n");
}

int main(int argc, char** argv) 
{	

   int ExitCode=EXIT_SUCCESS;
   while(1) {
      int option_index=0;
      int c=getopt_long(argc,argv,"hv",long_options,&option_index);
      if(c==-1) break;
      switch(c) {
      case 'h':
         Help(argv[0]);
         goto Exit1;
      case 'v':
         Version(argv[0]);
         goto Exit1;
      case '?':
         if(optopt=='?') {
            Help(argv[0]);
            goto Exit1;
         }
      default:
         ExitCode=EXIT_FAILURE;
         goto Exit1;
      }
   }

   if(argc==optind) {
      fprintf(stderr,"\nError, you must provide a file name and a database name as argument.\n");
      ExitCode=EXIT_FAILURE;
      goto Exit1;
       }
  
      char* filename= NULL;
      char* databasename= NULL;
 
   if(argc==(optind+1)) {
      filename=argv[optind];
      
   }

 else  if(argc==(optind+2)) {
       filename=argv[optind];
      databasename=argv[optind+1];
      
   }
   

   
   
   
   
   
   
   
   
// Open the file and starts the checking

   FILE* fp=fopen(filename,"rb");
   if(!fp) {
      fprintf(stderr,"\nError opening file `%s': %s.\n",filename,strerror(errno));
      ExitCode=EXIT_FAILURE;
      goto Exit1;
   }

   // Read magic bytes
   png_byte magic_png[8];
   if(fread(magic_png,1,sizeof(magic_png),fp)!=sizeof(magic_png)) {
      fprintf(stderr,"\nError, couldn't read %lu magic bytes from file %s.\n",(long unsigned)sizeof(magic_png),filename);
      ExitCode=EXIT_FAILURE;
      goto Exit2;
   }

   // Check for valid magic bytes
   if(!png_check_sig(magic_png,sizeof(magic_png))) {
      fprintf(stderr,"\nError, file `%s' is not a valid png image (invalid magic bytes).\n",filename);
      ExitCode=EXIT_FAILURE;
      goto Exit2;
   }

   // Create png read struct
   png_structp png_ptr=png_create_read_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
   if(!png_ptr) {
      fprintf(stderr,"\nError creating a png read struct.\n");
      ExitCode=EXIT_FAILURE;
      goto Exit2;
   }

   // Tell libpng to use fread with our file pointer
   png_init_io(png_ptr,fp);
   // Tell libpng we've already read the magic bytes
   png_set_sig_bytes(png_ptr,sizeof(magic_png));

   // Create png_info pointer
   png_infop info_ptr=png_create_info_struct(png_ptr);
   if(!info_ptr) {
      fprintf(stderr,"\nError, couldn't create an info pointer.\n");
      ExitCode=EXIT_FAILURE;
      goto Exit3;
   }

   // Read png info
   png_read_info(png_ptr,info_ptr);
   // DisplayImageInfo(png_ptr,info_ptr);

   // To check:
   // color_type==PNG_COLOR_TYPE_GRAY
   // channels==GRAY, PALETTE
   // bit_depth==8
   int width=png_get_image_width(png_ptr,info_ptr);
   int height=png_get_image_height(png_ptr,info_ptr);
   int bit_depth=png_get_bit_depth (png_ptr,info_ptr);
   int color_type=png_get_color_type (png_ptr,info_ptr);
   int channels=png_get_channels(png_ptr,info_ptr);
   unsigned rowbytes=png_get_rowbytes(png_ptr, info_ptr);
   if(bit_depth!=8) {
      fprintf(stderr,"\nError only images of bit depth 8 are supported. Here bit_depth=%d.\n",bit_depth);
      ExitCode=EXIT_FAILURE;
      goto Exit3;
   }
   if(color_type!=PNG_COLOR_TYPE_GRAY) {
      fprintf(stderr,"\nError only images with color type PNG_COLOR_TYPE_GRAY are supported.\n");
      ExitCode=EXIT_FAILURE;
      goto Exit3;
   }
   if(channels!=1) {
      fprintf(stderr,"\nError only images with GRAY channels (1) are supported. Here channels=%d\n",channels);
      ExitCode=EXIT_FAILURE;
      goto Exit3;
   }

   if(rowbytes!=width) {
      fprintf(stderr,"\nError, image rows need different number of bytes than its width: width=%d, nb of bytes=%d.\n",width,rowbytes);
      ExitCode=EXIT_FAILURE;
      goto Exit3;
   }

   // Allocate memory for image
   png_bytep png_row=(png_bytep)malloc(rowbytes*height);
   png_bytep* png_row_addresses=(png_bytep*)malloc(height*sizeof(png_bytep*));
   if(!png_row || !png_row_addresses) {
      fprintf(stderr,"\nError allocating %d bytes for `%s': %s.\n",rowbytes,filename,strerror(errno));
      ExitCode=EXIT_FAILURE;
      goto Exit3;
   }
   int i;
   for(i=0;i<height;i++) {
      png_row_addresses[i]=png_row+i*rowbytes;
   }
   // Read image in memory
   png_read_image(png_ptr,png_row_addresses);

   
  
   unsigned int x,y;
   long long unsigned GrayValue;
   long double J10, J01, J20, J11, J02, J30, J21, J12, J03;
   long long unsigned B00=0, B10=0, B01=0, B20=0, B11=0, B02=0, B30=0, B21=0, B12=0, B03=0;
   long long unsigned K10=0, K01=0, K20=0, K11=0, K02=0, K30=0, K21=0, K12=0, K03=0;
   //long long int R20=0, R11=0, R02=0, R30=0, R21=0, R03=0, R12=0;
   long double Mu20, Mu11, Mu02, Mu30, Mu21, Mu12, Mu03;
   
/* 
   First we compute B moments, which are a sum of discrete values 
   approximately, but not exactly, equal to J (raw) moments
*/

   for (x=0;x<width;++x)  {
      for (y=0;y<height;++y) {
         GrayValue=0xff-(long long unsigned)png_row[x*width+y];
         B00+=GrayValue;
         B10+=x*GrayValue;
         B01+=y*GrayValue;
         B20+=x*x*GrayValue;
         B11+=x*y*GrayValue;
         B02+=y*y*GrayValue;
         B30+=x*x*x*GrayValue;
         B21+=x*x*y*GrayValue;
         B12+=x*y*y*GrayValue;
         B03+=y*y*y*GrayValue;
      }

   }

   
// computing K moments (integers): J(p,q) = K(p,q) * some coefficient < 0
// K00 = B00 = J00

   K10=2*B10+B00;
   K01=2*B01+B00;
   K20=3*(B20+B10)+B00;
   K11=4*B11+2*(B10+B01)+B00;
   K02=3*(B02+B01)+B00;
   K30=4*(B30+B10)+6*B20+B00;
   K21=6*(B21+B11)+3*(B20+B10)+2*B01+B00;
   K12=6*(B12+B11)+3*(B02+B01)+2*B10+B00;
   K03=4*(B03+B01)+6*B02+B00;

   J10=(long double)K10/2.0;
   J01=(long double)K01/2.0;
   J20=(long double)K20/3.0;
   J11=(long double)K11/4.0;
   J02=(long double)K02/3.0;
   J30=(long double)K30/4.0;
   J21=(long double)K21/6.0;
   J12=(long double)K12/6.0;
   J03=(long double)K03/4.0;

  	// printf("K00 = %llu\n", B00);
	// printf("K10 = %llu\n", K10);
	// printf("K01 = %llu\n", K01);
	// printf("K20 = %llu\n", K20);
	// printf("K11 = %llu\n", K11);
	// printf("K02 = %llu\n", K02);
	// printf("K30 = %llu\n", K30);
	// printf("K21 = %llu\n", K21);
	// printf("K12 = %llu\n", K12); 
	// printf("K03 = %llu\n", K03); 


  	// printf("J00 = %.20llu\n", B00);
	// printf("J10 = %.20Lf\n", J10);
	// printf("J01 = %.20Lf\n", J01);
	// printf("J20 = %.20Lf\n", J20);
	// printf("J11 = %.20Lf\n", J11);
	// printf("J02 = %.20Lf\n", J02);
	// printf("J30 = %.20Lf\n", J30);
	// printf("J21 = %.20Lf\n", J21);
	// printf("J12 = %.20Lf\n", J12); 
	// printf("J03 = %.20Lf\n", J03); 


// computing R moments (integers): µ(p,q) = R(p,q) * some coefficient < 0
// R00 = µ00 = J00 = B00
// R10 = R01 = µ10 = µ01 = 0

   // R20=4*K20*B00-3*K10*K10;
   // R11=K11*B00-K10*K01;
   // R02=4*K02*B00-3*K01*K01;
   // R30=K10*(K10*K10-K20*B00)+B00*(K30*B00-K10*K20);
   // R21=3*K10*(K10*K01-K11*B00)+2*B00*(K21*B00-K20*K01);
   // R12=3*K01*(K10*K01-K11*B00)+2*B00*(K12*B00-K02*K10);
   // R03=K01*(K01*K01-K02*B00)+B00*(K03*B00-K01*K02);

   // Mu20=(long double)R20/(12*B00);
   // Mu11=(long double)R11/(4*B00);
   // Mu02=(long double)R02/(12*B00);
   // Mu30=(long double)R30/(4*B00*B00);
   // Mu21=(long double)R21/(12*B00*B00);
   // Mu12=(long double)R12/(12*B00*B00);
   // Mu03=(long double)R03/(4*B00*B00);
   Mu20=(J20*B00-J10*J10)/B00;
   Mu02=(J02*B00-J01*J01)/B00;
   Mu11=(J11*B00-J01*J10)/B00;
   Mu30=(J30*B00*B00-3*J10*J20*B00+2*J10*J10*J10)/(B00*B00);
   Mu03=(J03*B00*B00-3*J01*J02*B00+2*J01*J01*J01)/(B00*B00);
   Mu21=(2*J10*J10*J01-J01*J20*B00-2*J10*J11*B00+J21*B00*B00)/(B00*B00);
   Mu12=(2*J01*J01*J10-J10*J02*B00-2*J01*J11*B00+J12*B00*B00)/(B00*B00);

	// printf("Mu20 = %.10Lf\n", Mu20);
	// printf("Mu11 = %.10Lf\n", Mu11);
	// printf("Mu02 = %.10Lf\n", Mu02);
	// printf("Mu30 = %.10Lf\n", Mu30);
	// printf("Mu21 = %.10Lf\n", Mu21);
	// printf("Mu12 = %.10Lf\n", Mu12); 
	// printf("Mu03 = %.10Lf\n", Mu03); 
   
// And now we compute the M moments, which are the N moments * 255 ^ (p+q)/2 
// With N(p,q) = µ(p,q) / [ µ(0,0) ^ (1 + (p+q)/2) ] 
      
   long double M20, M11, M02, M30, M21, M12, M03;
   long double normcoeff=pow(255,1.5);
   long double normB00=pow(B00,2.5);

   // M20 = 255* (1.0/(12.0*B00)) * R20 / (B00*B00);
   // M11 = 255* (1.0/(4.0*B00))  * R11 / (B00*B00);
   // M02 = 255* (1.0/(12.0*B00)) * R02 / (B00*B00);
   // M30 = normcoeff * (1.0/(4.0 *B00*B00)) * R30 / pow(B00,2.5);
   // M21 = normcoeff * (1.0/(12.0*B00*B00)) * R21 / pow(B00,2.5);
   // M12 = normcoeff * (1.0/(12.0*B00*B00)) * R12 / pow(B00,2.5);
   // M03 = normcoeff * (1.0/(4.0 *B00*B00)) * R03 / pow(B00,2.5);
   M20=255*Mu20/(B00*B00);
   M11=255*Mu11/(B00*B00);
   M02=255*Mu02/(B00*B00);
   M30=normcoeff*Mu30/normB00;
   M03=normcoeff*Mu03/normB00;
   M21=normcoeff*Mu21/normB00;
   M12=normcoeff*Mu12/normB00;


// And at last, the seven invariant moments of the analysed image

	long double I_UnknownFig[7]={0}, I_ReferenceFig[7]={0};


	I_UnknownFig[0]=M20+M02;

	I_UnknownFig[1]=(M20-M02)*(M20-M02)+(2*M11)*(2*M11);

	I_UnknownFig[2]=(M30-3*M12)*(M30-3*M12)+(3*M21-M03)*(3*M21-M03);

	I_UnknownFig[3]=(M30+M12)*(M30+M12)+(M21+M03)*(M21+M03);

	I_UnknownFig[4]= (M30-3*M12)*(M30+M12)*( (M30+M12)*(M30+M12) - 3*(M21+M03)*(M21+M03) ) + (3*M21-M03)*(M21+M03)*( 3*(M30+M12)*(M30+M12) - (M21+M03)*(M21+M03) );

	I_UnknownFig[5]= (M20-M02)*( (M30+M12)*(M30+M12) - (M21+M03)*(M21+M03) ) + 4*M11*(M30+M12)*(M21+M03);

	I_UnknownFig[6]= (3*M21-M03)*(M30+M12)*( (M30+M12)*(M30+M12) - 3*(M21+M03)*(M21+M03) ) + (M30-3*M12)*(M21+M03)*( 3*(M30+M12)*(M30+M12) - (M21+M03)*(M21+M03) );

      
      x=0;
      printf("\n\n");
      do {
         printf("I%d = %.25Lf\n",x+1,I_UnknownFig[x]);
         ++x;
      } while(x<7);
      printf("\n");
      
     
     
     
     
if (databasename != NULL){
     

// In this part we determine which figure is the closest to the input figure. 


int k=0;
long double DistanceToRef=0.0, DistanceMin=0.0;
char ReferenceFig[50]={0}, ClosestFig[50]={0};

FILE* database = fopen(databasename, "r+");

while (fgetc(database) != EOF)    // as long as there is a character to read in the database do : 
	{	
		fseek(database,-1,SEEK_CUR);  // go at the first character position off the current line
		
		fscanf(database, "%s %Lf %Lf %Lf %Lf %Lf %Lf %Lf\n", ReferenceFig, &I_ReferenceFig[0], &I_ReferenceFig[1], &I_ReferenceFig[2], &I_ReferenceFig[3], &I_ReferenceFig[4], &I_ReferenceFig[5], &I_ReferenceFig[6]);            // read the name of the reference figure and store its 7 I moments in I_Reference[] 

		DistanceToRef=0.0;
		for (i=0;i<7;++i)
			{
				DistanceToRef += (I_UnknownFig[i]-I_ReferenceFig[i])*(I_UnknownFig[i]-I_ReferenceFig[i]);  // compute the Distance of the I of your image compared to the current reference image. 
			}
		if(k==0) {
				DistanceMin = DistanceToRef;
				strcpy(ClosestFig,ReferenceFig);   // if this is the first comparison (using the first line) of the database, the closest figure we know of your image must be this reference image 
				k=1;
			}	
		else if (DistanceToRef < DistanceMin)
			{
				DistanceMin=DistanceToRef;
				strcpy(ClosestFig,ReferenceFig);
			}
		printf("\nDistance to %s : %.25Lf", ReferenceFig, DistanceToRef);
	}
		
printf("\n\nThis must be a %s!\n\nDo you think this is something else? [Y/N]     ",ClosestFig);
char answer[50];
scanf("%s",answer);                                               // in this part we store the moments of the input figure if the user decides so and names this new figure.
if (answer[0]=='Y' || answer[0]=='y')
{		printf("What is it then?   ");
		char newfigurename[50];
		scanf("%s",newfigurename);
		fseek(database, 0, SEEK_END);
		fprintf(database, "%s %.25Lf %.25Lf %.25Lf %.25Lf %.25Lf %.25Lf %.25Lf\n", newfigurename, I_UnknownFig[0], I_UnknownFig[1], I_UnknownFig[2], I_UnknownFig[3], I_UnknownFig[4], I_UnknownFig[5], I_UnknownFig[6]);
		fclose(database);
}
else {fclose(database);}
     
     
}
 
//  And then there are those last lines, containing the exits for the goto placed in the initial checkings
 
   free(png_row);
Exit3:
   png_destroy_read_struct(&png_ptr,&info_ptr,NULL);
Exit2:
   fclose(fp);
Exit1:
   exit(ExitCode);



}













