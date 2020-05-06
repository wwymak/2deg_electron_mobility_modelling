#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define YES 1
#define NO 0


int getline_from_file(FILE *, char *,char *);
int remove_comments_from_linestring(char *); /* comments begin with a #, but \# means that the # is required */
int getsection_from_file(FILE *, char *, char *, int *, char *);
int get_number_of_numbers_in_a_string(char *);
int search_string_and_update_value(char *, char *, char *);




int getline_from_file (FILE *stream, char *destination_string,  char *ptr_to_last_char)
{  /* return value is the last character of the line.. carriage return or EOF */

 char ch;
 int found_first_nonwhite_char=NO, i,j=0;


  do { ch=fgetc(stream); /* fscanf will not read EOF; need to read the return value of the fscanf for that */
       switch(ch)
	       { case ' ' : /* white characters : ignore if they occur before any non white char */
		     case '\t': if(found_first_nonwhite_char) destination_string[j++]=ch;
			            break;

			 case '\n': /* these are not to be copied to datastring  */
			 case EOF : *ptr_to_last_char =ch;
			            destination_string [j]='\0';
		   	            break;

		     default  : found_first_nonwhite_char=YES;
			            destination_string[j++]=ch;
			            break;
		    }

     }while( (ch!=EOF) && (ch!='\n') );

 return strlen(destination_string);
 }



int remove_comments_from_linestring(char *linestring)

{  /* comments begin with a #, but \# means that the # is required */


  int i,j,n;
  n=strlen(linestring);

  j=0;
  for(i=1; i<=n; i++)  /* the first char cannot be #, then this prog would not have been called anyway */
      {                /* ensure that '\0' is copied correctly */

          if(linestring[i]=='#')
			{

			  switch(linestring[i-1])
			         {
					   case '\\':linestring[j]=linestring[i];
					             break;

				       default :linestring[i]='\0';
					            return j++; /* j is the number of characters retained in the string */
				                break;
			         }

		    } else{ linestring[++j]=linestring[i];}


	  }

return 0;

}


 int getsection_from_file(FILE *stream, char *section_name, char *destination_string, int *errcode_ptr, char *errmsg)

 {
  int i,this_is_comment_line=YES,j1,j2,j3,linenumber,section_complete, found_begin_section, found_end_section;
  char inputline[65536], ch, *pch;

  /* scan the FILE for \begin{section_name} */
   rewind(stream);


   section_complete=NO;
   found_begin_section=NO;
   found_end_section=NO;
   linenumber=0;

   printf("looking for : %s\n", section_name);

  do{



    if(!found_begin_section)
	  {
     i=getline_from_file (stream,inputline,&ch);
	 remove_comments_from_linestring(inputline); linenumber++;
     /* printf("getline_from_file  returned %s after line %d\n",(ch==EOF)? "EOF":"newline",linenumber); */

      if((i>0) && (inputline[0]!='#')) /* it is not a blank line  or a comment line */
	    {
		  pch = strstr(inputline,"\\begin");  /* some section is starting */
	      if(pch){

			       pch=strstr(pch+5,"{");
			       if(pch){
					       pch=strstr(pch+1,section_name);
					       if(pch){
							       printf("line %d : %s begins\n",linenumber,section_name);
							       found_begin_section=YES;

						          }
	  				       }
			      }


		 }


       }


	   if(found_begin_section){

                                i=getline_from_file (stream,inputline,&ch); linenumber++;
                                remove_comments_from_linestring(inputline);
								if(i) /* it is not a blank line */
	                               {pch = strstr(inputline,"\\end");  /* keep looking for \end{section} */
	                                if(pch){ pch=strstr(pch,"{");
			                                if(pch){pch=strstr(pch,section_name);
					                                if(pch)
													  { printf("line %d : %s ends\n",linenumber,section_name);
													   found_end_section=YES;
					                                  }
													}


				                            }


                                   }



                                 /* concatenate the input to destination_string */
	                             if(!found_end_section){ /* check this is not a comment line */
									                     /* whose first nonwhite char is # */
														i=0;
														while( (inputline[i] == '\t') || (inputline[i] == ' ') ){i++;}

														if(inputline[i] != '#'){
								                                                strcat(destination_string,inputline);
														                        strcat(destination_string," ");
																			   }

													   }


						       }





      section_complete = ((found_begin_section==YES) && (found_end_section==YES))? YES : NO;


     }while( (ch!=EOF) && (section_complete==NO) );




      if(section_complete==NO){*errcode_ptr=-1; sprintf(errmsg,"section not found or complete");}


 rewind(stream);
 if(section_complete==YES) {return strlen(destination_string);
                           }else {return -1;
						         }


 }




 int search_string_and_update_value(char *varname,  char *inputstring, char *returnvaluestring)
{
 /* return 1 : if a value was found
    return 0 : no match was found for the varname string
 */


 int n=0,i,j;
 char *pch, *localcopy;

 /* make a local copy so that the original string is untouched */

 i=strlen(inputstring);
 localcopy=(char *) malloc((i+1)*sizeof(char));
 strcpy(localcopy,inputstring);

  pch = strtok (localcopy," =\t");  /* the delimitor can be space or = or \t */
  while (pch != NULL)
  {
    i=strcmp(pch,varname);           /* if a token matches varname, then next token must be its value */
    pch = strtok (NULL, " =\t");
	if((i==0) && (strlen(pch) > 0))
	  {strcpy(returnvaluestring,pch);
	   return 1;
	  }
   }


 free(localcopy);

 return 0;

}





int get_number_of_numbers_in_a_string(char *inputstring)
{

 int n=0,i,j;
 char *pch, *localcopy;

 /* make a local copy so that the original string is untouched */

 i=strlen(inputstring);
 localcopy=(char *) malloc((i+1)*sizeof(char));
 strcpy(localcopy,inputstring);

  pch = strtok (localcopy," ,;\t");  /* the delimitor can be space or , or ; or \t */
  while (pch != NULL)
  {
    n++;
    pch = strtok (NULL, " ,;\t");
  }


free(localcopy);
return n;

}

