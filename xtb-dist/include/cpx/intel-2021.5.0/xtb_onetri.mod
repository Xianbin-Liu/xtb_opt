  +@  �   k820309    ,          2021.5.0    ֣�f                                                                                                          
       ../src/onetri.f90 XTB_ONETRI              ONETRI                                                     
                                                           
       WP                                                     
       BLAS_SYMM BLAS_GEMM               @�                                       u #SGEMM    #DGEMM    #CGEMM     #ZGEMM .   #SCGEMM <   #DZGEMM J   #         @     @                                	               #TRANSA    #TRANSB    #M    #N    #K 	   #ALPHA 
   #A    #LDA    #B    #LDB    #BETA    #C    #LDC              
                                                                     
                                                                     
                                                      
                                                      
                                 	                     
                                
     	             B  
                                                    	      p        5 O p        p          5 O p          1     5 O p          1                             
                                                  B  
                                                    	      p        5 O p 
       p          5 O p 
         1     5 O p 
         1                             
                                                     
                                     	             B  
                                                   	       p        5 O p        p          5 O p          1     5 O p          1                             
                                           #         @     @                                	               #TRANSA    #TRANSB    #M    #N    #K    #ALPHA    #A    #LDA    #B    #LDB    #BETA    #C    #LDC              
                                                                     
                                                                     
                                                      
                                                      
                                                      
                                     
             B  
                                                    
      p        5 O p        p          5 O p          1     5 O p          1                             
                                                  B  
                                                    
      p        5 O p 
       p          5 O p 
         1     5 O p 
         1                             
                                                     
                                     
             B  
                                                   
       p        5 O p        p          5 O p          1     5 O p          1                             
                                           #         @     @                                 	               #TRANSA !   #TRANSB "   #M #   #N $   #K %   #ALPHA &   #A '   #LDA (   #B )   #LDB *   #BETA +   #C ,   #LDC -             
                                !                                     
                                "                                     
                                 #                     
                                 $                     
                                 %                     
                                &                  B  
                                '                          p        5 O p        p          5 O p          1     5 O p          1                             
                                (                  B  
                                )                          p        5 O p 
       p          5 O p 
         1     5 O p 
         1                             
                                *                     
                                +                  B  
                               ,                     	      p        5 O p        p          5 O p          1     5 O p          1                             
                                -           #         @     @                           .     	               #TRANSA /   #TRANSB 0   #M 1   #N 2   #K 3   #ALPHA 4   #A 5   #LDA 6   #B 7   #LDB 8   #BETA 9   #C :   #LDC ;             
                                /                                     
                                0                                     
                                 1                     
                                 2                     
                                 3                     
                                4                  B  
                                5                     
     p        5 O p        p          5 O p          1     5 O p          1                             
                                6                  B  
                                7                          p        5 O p 
       p          5 O p 
         1     5 O p 
         1                             
                                8                     
                                9                  B  
                               :                           p        5 O p        p          5 O p          1     5 O p          1                             
                                ;           #         @     @                           <     	               #TRANSA =   #TRANSB >   #M ?   #N @   #K A   #ALPHA B   #A C   #LDA D   #B E   #LDB F   #BETA G   #C H   #LDC I             
                                =                                     
                                >                                     
                                 ?                     
                                 @                     
                                 A                     
                                B                  B  
                                C                    	      p        5 O p        p          5 O p          1     5 O p          1                             
                                D                  B  
                                E                          p        5 O p 
       p          5 O p 
         1     5 O p 
         1                             
                                F                     
                                G                  B  
                               H                           p        5 O p        p          5 O p          1     5 O p          1                             
                                I           #         @     @                           J     	               #TRANSA K   #TRANSB L   #M M   #N N   #K O   #ALPHA P   #A Q   #LDA R   #B S   #LDB T   #BETA U   #C V   #LDC W             
                                K                                     
                                L                                     
                                 M                     
                                 N                     
                                 O                     
                                P                  B  
                                Q                    
      p        5 O p        p          5 O p          1     5 O p          1                             
                                R                  B  
                                S                          p        5 O p 
       p          5 O p 
         1     5 O p 
         1                             
                                T                     
                                U                  B  
                               V                           p        5 O p        p          5 O p          1     5 O p          1                             
                                W                         @�                                       u #SSYMM X   #DSYMM e   #CSYMM r   #ZSYMM    #         @     @                           X     	               #SIDE Y   #UPLO Z   #M [   #N \   #ALPHA ]   #A ^   #LDA _   #B `   #LDB a   #BETA b   #C c   #LDC d             
                                Y                                     
                                Z                                     
                                 [                     
                                 \                     
                                ]     	             B  
                                ^                    	 #     p        5 O p        p          5 O p          1     5 O p          1                             
                                _                  B  
                                `                    	 $     p        5 O p 	       p          5 O p 	         1     5 O p 	         1                             
                                a                     
                                b     	             B  
                               c                    	 %      p        5 O p        p          5 O p          1     5 O p          1                             
                                d           #         @     @                           e     	               #SIDE f   #UPLO g   #M h   #N i   #ALPHA j   #A k   #LDA l   #B m   #LDB n   #BETA o   #C p   #LDC q             
                                f                                     
                                g                                     
                                 h                     
                                 i                     
                                j     
             B  
                                k                    
 &     p        5 O p        p          5 O p          1     5 O p          1                             
                                l                  B  
                                m                    
 '     p        5 O p 	       p          5 O p 	         1     5 O p 	         1                             
                                n                     
                                o     
             B  
                               p                    
 (      p        5 O p        p          5 O p          1     5 O p          1                             
                                q           #         @     @                           r     	               #SIDE s   #UPLO t   #M u   #N v   #ALPHA w   #A x   #LDA y   #B z   #LDB {   #BETA |   #C }   #LDC ~             
                                s                                     
                                t                                     
                                 u                     
                                 v                     
                                w                  B  
                                x                     )     p        5 O p        p          5 O p          1     5 O p          1                             
                                y                  B  
                                z                     *     p        5 O p 	       p          5 O p 	         1     5 O p 	         1                             
                                {                     
                                |                  B  
                               }                     +      p        5 O p        p          5 O p          1     5 O p          1                             
                                ~           #         @     @                                	               #SIDE �   #UPLO �   #M �   #N �   #ALPHA �   #A �   #LDA �   #B �   #LDB �   #BETA �   #C �   #LDC �             
                                �                                     
                                �                                     
                                 �                     
                                 �                     
                                �                  B  
                                �                     ,     p        5 O p        p          5 O p          1     5 O p          1                             
                                �                  B  
                                �                     -     p        5 O p 	       p          5 O p 	         1     5 O p 	         1                             
                                �                     
                                �                  B  
                               �                     .      p        5 O p        p          5 O p          1     5 O p          1                             
                                �           #         @                                   �                    #ITY �   #S �   #S1 �   #ARRAY �   #N �   #IVAL �             
  @                               �                  @  
  @                              �                    
    p          1     1                          @  
D @                              �                    
     p          1     1                            
@ @                              �                    
      p        5 � p        r �   p          5 � p        r �     5 � p        r �       5 � p        r �     5 � p        r �                               
@ @                               �                     
@ @                               �              �   %      fn#fn     �      b   uapp(XTB_ONETRI    �   @   J  XTB_BLOWSY "     C   J  XTB_MCTC_ACCURACY    _  T   J  XTB_MCTC_BLAS ,   �  �       gen@BLAS_GEMM+XTB_MCTC_BLAS +   7  �      SGEMM+XTB_MCTC_BLAS_LEVEL3 2   �  P   a   SGEMM%TRANSA+XTB_MCTC_BLAS_LEVEL3 2   A  P   a   SGEMM%TRANSB+XTB_MCTC_BLAS_LEVEL3 -   �  @   a   SGEMM%M+XTB_MCTC_BLAS_LEVEL3 -   �  @   a   SGEMM%N+XTB_MCTC_BLAS_LEVEL3 -     @   a   SGEMM%K+XTB_MCTC_BLAS_LEVEL3 1   Q  @   a   SGEMM%ALPHA+XTB_MCTC_BLAS_LEVEL3 -   �  �   a   SGEMM%A+XTB_MCTC_BLAS_LEVEL3 /   m  @   a   SGEMM%LDA+XTB_MCTC_BLAS_LEVEL3 -   �  �   a   SGEMM%B+XTB_MCTC_BLAS_LEVEL3 /   �  @   a   SGEMM%LDB+XTB_MCTC_BLAS_LEVEL3 0   �  @   a   SGEMM%BETA+XTB_MCTC_BLAS_LEVEL3 -   	  �   a   SGEMM%C+XTB_MCTC_BLAS_LEVEL3 /   �  @   a   SGEMM%LDC+XTB_MCTC_BLAS_LEVEL3 +   %  �      DGEMM+XTB_MCTC_BLAS_LEVEL3 2   �  P   a   DGEMM%TRANSA+XTB_MCTC_BLAS_LEVEL3 2   /	  P   a   DGEMM%TRANSB+XTB_MCTC_BLAS_LEVEL3 -   	  @   a   DGEMM%M+XTB_MCTC_BLAS_LEVEL3 -   �	  @   a   DGEMM%N+XTB_MCTC_BLAS_LEVEL3 -   �	  @   a   DGEMM%K+XTB_MCTC_BLAS_LEVEL3 1   ?
  @   a   DGEMM%ALPHA+XTB_MCTC_BLAS_LEVEL3 -   
  �   a   DGEMM%A+XTB_MCTC_BLAS_LEVEL3 /   [  @   a   DGEMM%LDA+XTB_MCTC_BLAS_LEVEL3 -   �  �   a   DGEMM%B+XTB_MCTC_BLAS_LEVEL3 /   w  @   a   DGEMM%LDB+XTB_MCTC_BLAS_LEVEL3 0   �  @   a   DGEMM%BETA+XTB_MCTC_BLAS_LEVEL3 -   �  �   a   DGEMM%C+XTB_MCTC_BLAS_LEVEL3 /   �  @   a   DGEMM%LDC+XTB_MCTC_BLAS_LEVEL3 +     �      CGEMM+XTB_MCTC_BLAS_LEVEL3 2   �  P   a   CGEMM%TRANSA+XTB_MCTC_BLAS_LEVEL3 2     P   a   CGEMM%TRANSB+XTB_MCTC_BLAS_LEVEL3 -   m  @   a   CGEMM%M+XTB_MCTC_BLAS_LEVEL3 -   �  @   a   CGEMM%N+XTB_MCTC_BLAS_LEVEL3 -   �  @   a   CGEMM%K+XTB_MCTC_BLAS_LEVEL3 1   -  @   a   CGEMM%ALPHA+XTB_MCTC_BLAS_LEVEL3 -   m  �   a   CGEMM%A+XTB_MCTC_BLAS_LEVEL3 /   I  @   a   CGEMM%LDA+XTB_MCTC_BLAS_LEVEL3 -   �  �   a   CGEMM%B+XTB_MCTC_BLAS_LEVEL3 /   e  @   a   CGEMM%LDB+XTB_MCTC_BLAS_LEVEL3 0   �  @   a   CGEMM%BETA+XTB_MCTC_BLAS_LEVEL3 -   �  �   a   CGEMM%C+XTB_MCTC_BLAS_LEVEL3 /   �  @   a   CGEMM%LDC+XTB_MCTC_BLAS_LEVEL3 +     �      ZGEMM+XTB_MCTC_BLAS_LEVEL3 2   �  P   a   ZGEMM%TRANSA+XTB_MCTC_BLAS_LEVEL3 2     P   a   ZGEMM%TRANSB+XTB_MCTC_BLAS_LEVEL3 -   [  @   a   ZGEMM%M+XTB_MCTC_BLAS_LEVEL3 -   �  @   a   ZGEMM%N+XTB_MCTC_BLAS_LEVEL3 -   �  @   a   ZGEMM%K+XTB_MCTC_BLAS_LEVEL3 1     @   a   ZGEMM%ALPHA+XTB_MCTC_BLAS_LEVEL3 -   [  �   a   ZGEMM%A+XTB_MCTC_BLAS_LEVEL3 /   7  @   a   ZGEMM%LDA+XTB_MCTC_BLAS_LEVEL3 -   w  �   a   ZGEMM%B+XTB_MCTC_BLAS_LEVEL3 /   S  @   a   ZGEMM%LDB+XTB_MCTC_BLAS_LEVEL3 0   �  @   a   ZGEMM%BETA+XTB_MCTC_BLAS_LEVEL3 -   �  �   a   ZGEMM%C+XTB_MCTC_BLAS_LEVEL3 /   �  @   a   ZGEMM%LDC+XTB_MCTC_BLAS_LEVEL3 ,   �  �      SCGEMM+XTB_MCTC_BLAS_LEVEL3 3   �  P   a   SCGEMM%TRANSA+XTB_MCTC_BLAS_LEVEL3 3   �  P   a   SCGEMM%TRANSB+XTB_MCTC_BLAS_LEVEL3 .   I  @   a   SCGEMM%M+XTB_MCTC_BLAS_LEVEL3 .   �  @   a   SCGEMM%N+XTB_MCTC_BLAS_LEVEL3 .   �  @   a   SCGEMM%K+XTB_MCTC_BLAS_LEVEL3 2   	  @   a   SCGEMM%ALPHA+XTB_MCTC_BLAS_LEVEL3 .   I  �   a   SCGEMM%A+XTB_MCTC_BLAS_LEVEL3 0   %  @   a   SCGEMM%LDA+XTB_MCTC_BLAS_LEVEL3 .   e  �   a   SCGEMM%B+XTB_MCTC_BLAS_LEVEL3 0   A  @   a   SCGEMM%LDB+XTB_MCTC_BLAS_LEVEL3 1   �  @   a   SCGEMM%BETA+XTB_MCTC_BLAS_LEVEL3 .   �  �   a   SCGEMM%C+XTB_MCTC_BLAS_LEVEL3 0   �  @   a   SCGEMM%LDC+XTB_MCTC_BLAS_LEVEL3 ,   �  �      DZGEMM+XTB_MCTC_BLAS_LEVEL3 3   �   P   a   DZGEMM%TRANSA+XTB_MCTC_BLAS_LEVEL3 3   �   P   a   DZGEMM%TRANSB+XTB_MCTC_BLAS_LEVEL3 .   7!  @   a   DZGEMM%M+XTB_MCTC_BLAS_LEVEL3 .   w!  @   a   DZGEMM%N+XTB_MCTC_BLAS_LEVEL3 .   �!  @   a   DZGEMM%K+XTB_MCTC_BLAS_LEVEL3 2   �!  @   a   DZGEMM%ALPHA+XTB_MCTC_BLAS_LEVEL3 .   7"  �   a   DZGEMM%A+XTB_MCTC_BLAS_LEVEL3 0   #  @   a   DZGEMM%LDA+XTB_MCTC_BLAS_LEVEL3 .   S#  �   a   DZGEMM%B+XTB_MCTC_BLAS_LEVEL3 0   /$  @   a   DZGEMM%LDB+XTB_MCTC_BLAS_LEVEL3 1   o$  @   a   DZGEMM%BETA+XTB_MCTC_BLAS_LEVEL3 .   �$  �   a   DZGEMM%C+XTB_MCTC_BLAS_LEVEL3 0   �%  @   a   DZGEMM%LDC+XTB_MCTC_BLAS_LEVEL3 ,   �%  l       gen@BLAS_SYMM+XTB_MCTC_BLAS +   7&  �      SSYMM+XTB_MCTC_BLAS_LEVEL3 0   �&  P   a   SSYMM%SIDE+XTB_MCTC_BLAS_LEVEL3 0   6'  P   a   SSYMM%UPLO+XTB_MCTC_BLAS_LEVEL3 -   �'  @   a   SSYMM%M+XTB_MCTC_BLAS_LEVEL3 -   �'  @   a   SSYMM%N+XTB_MCTC_BLAS_LEVEL3 1   (  @   a   SSYMM%ALPHA+XTB_MCTC_BLAS_LEVEL3 -   F(  �   a   SSYMM%A+XTB_MCTC_BLAS_LEVEL3 /   ")  @   a   SSYMM%LDA+XTB_MCTC_BLAS_LEVEL3 -   b)  �   a   SSYMM%B+XTB_MCTC_BLAS_LEVEL3 /   >*  @   a   SSYMM%LDB+XTB_MCTC_BLAS_LEVEL3 0   ~*  @   a   SSYMM%BETA+XTB_MCTC_BLAS_LEVEL3 -   �*  �   a   SSYMM%C+XTB_MCTC_BLAS_LEVEL3 /   �+  @   a   SSYMM%LDC+XTB_MCTC_BLAS_LEVEL3 +   �+  �      DSYMM+XTB_MCTC_BLAS_LEVEL3 0   �,  P   a   DSYMM%SIDE+XTB_MCTC_BLAS_LEVEL3 0   �,  P   a   DSYMM%UPLO+XTB_MCTC_BLAS_LEVEL3 -   )-  @   a   DSYMM%M+XTB_MCTC_BLAS_LEVEL3 -   i-  @   a   DSYMM%N+XTB_MCTC_BLAS_LEVEL3 1   �-  @   a   DSYMM%ALPHA+XTB_MCTC_BLAS_LEVEL3 -   �-  �   a   DSYMM%A+XTB_MCTC_BLAS_LEVEL3 /   �.  @   a   DSYMM%LDA+XTB_MCTC_BLAS_LEVEL3 -   /  �   a   DSYMM%B+XTB_MCTC_BLAS_LEVEL3 /   �/  @   a   DSYMM%LDB+XTB_MCTC_BLAS_LEVEL3 0   !0  @   a   DSYMM%BETA+XTB_MCTC_BLAS_LEVEL3 -   a0  �   a   DSYMM%C+XTB_MCTC_BLAS_LEVEL3 /   =1  @   a   DSYMM%LDC+XTB_MCTC_BLAS_LEVEL3 +   }1  �      CSYMM+XTB_MCTC_BLAS_LEVEL3 0   ,2  P   a   CSYMM%SIDE+XTB_MCTC_BLAS_LEVEL3 0   |2  P   a   CSYMM%UPLO+XTB_MCTC_BLAS_LEVEL3 -   �2  @   a   CSYMM%M+XTB_MCTC_BLAS_LEVEL3 -   3  @   a   CSYMM%N+XTB_MCTC_BLAS_LEVEL3 1   L3  @   a   CSYMM%ALPHA+XTB_MCTC_BLAS_LEVEL3 -   �3  �   a   CSYMM%A+XTB_MCTC_BLAS_LEVEL3 /   h4  @   a   CSYMM%LDA+XTB_MCTC_BLAS_LEVEL3 -   �4  �   a   CSYMM%B+XTB_MCTC_BLAS_LEVEL3 /   �5  @   a   CSYMM%LDB+XTB_MCTC_BLAS_LEVEL3 0   �5  @   a   CSYMM%BETA+XTB_MCTC_BLAS_LEVEL3 -   6  �   a   CSYMM%C+XTB_MCTC_BLAS_LEVEL3 /   �6  @   a   CSYMM%LDC+XTB_MCTC_BLAS_LEVEL3 +    7  �      ZSYMM+XTB_MCTC_BLAS_LEVEL3 0   �7  P   a   ZSYMM%SIDE+XTB_MCTC_BLAS_LEVEL3 0   8  P   a   ZSYMM%UPLO+XTB_MCTC_BLAS_LEVEL3 -   o8  @   a   ZSYMM%M+XTB_MCTC_BLAS_LEVEL3 -   �8  @   a   ZSYMM%N+XTB_MCTC_BLAS_LEVEL3 1   �8  @   a   ZSYMM%ALPHA+XTB_MCTC_BLAS_LEVEL3 -   /9  �   a   ZSYMM%A+XTB_MCTC_BLAS_LEVEL3 /   :  @   a   ZSYMM%LDA+XTB_MCTC_BLAS_LEVEL3 -   K:  �   a   ZSYMM%B+XTB_MCTC_BLAS_LEVEL3 /   ';  @   a   ZSYMM%LDB+XTB_MCTC_BLAS_LEVEL3 0   g;  @   a   ZSYMM%BETA+XTB_MCTC_BLAS_LEVEL3 -   �;  �   a   ZSYMM%C+XTB_MCTC_BLAS_LEVEL3 /   �<  @   a   ZSYMM%LDC+XTB_MCTC_BLAS_LEVEL3    �<  |       ONETRI    ?=  @   a   ONETRI%ITY    =  �   a   ONETRI%S    >  �   a   ONETRI%S1    �>  $  a   ONETRI%ARRAY    �?  @   a   ONETRI%N    �?  @   a   ONETRI%IVAL 