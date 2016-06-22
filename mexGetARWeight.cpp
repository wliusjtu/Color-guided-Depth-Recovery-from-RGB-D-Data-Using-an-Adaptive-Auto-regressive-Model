#include<mex.h>
#include<math.h>

void GetDepthBlock(int cur_i, int cur_j, int m, int n, int winR, double *Image, double *Block, double *Coor, int InvFlag);
void GetColorPatch(int cur_i, int cur_j, int mExtend, int n, int winR, int rPatch, double *Image, double *Block, double *Coor, int InvFlag);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *CR, *CG, *CB, *DepthGuide;
    double sigma1, sigma2, sigma4;
    int m, n, rWin, rPatch, winSize, blockSize;
    
    CR = mxGetPr(prhs[0]); // the input color image MUST be padarray 2*rPatch pixels in each dimension
    CG = mxGetPr(prhs[1]);
    CB = mxGetPr(prhs[2]);
    DepthGuide = mxGetPr(prhs[3]); // the guide depth map
    rWin = (int)mxGetScalar(prhs[4]); // the radius of the AR neighborhood system N(x) in Eq.(5)
    rPatch = (int)mxGetScalar(prhs[5]); // the radius of the patch to compute the shape-based color term of the AR coefficient.
    sigma1 = mxGetScalar(prhs[6]);
    sigma2 = mxGetScalar(prhs[7]);
    sigma4 = mxGetScalar(prhs[8]);
  
    
    m = mxGetM(prhs[3]); // row number
    n = mxGetN(prhs[3]); // column number
    int mExtend = m + 2*rPatch;
    
    
    winSize = 2*rWin + 1; //window size of the AR window
    blockSize = 2*rPatch + 1; // block size to compute the shape-based color weight
    int winSizeExtend = winSize + 2*rPatch;
    
    /************* the output Q matrix **********************/
    plhs[0] = mxCreateSparse(m*n, m*n, winSize*winSize*m*n, mxREAL); 
    double *ARWeightPr = mxGetPr(plhs[0]); // the output depth weight sparse matrix for data term
    mwIndex *ARWeightIr = mxGetIr(plhs[0]);
    mwIndex *ARWeightJc = mxGetJc(plhs[0]);

    
    /************ compute weight *************/
    int i, j, s, t, p, q, counter;
    int InvFlag = -1; // if the pixle in the window is outside the image, set the corresponding coordinate in Coor as InvFlag.
    double diffR, diffG, diffB, diffColor, diffDepth, WeightSum;
    double *BlockR = mxGetPr(mxCreateDoubleMatrix(winSizeExtend, winSizeExtend, mxREAL));
    double *BlockG = mxGetPr(mxCreateDoubleMatrix(winSizeExtend, winSizeExtend, mxREAL));
    double *BlockB = mxGetPr(mxCreateDoubleMatrix(winSizeExtend, winSizeExtend, mxREAL));
    double *BlockDepth = mxGetPr(mxCreateDoubleMatrix(winSize, winSize, mxREAL));
    double *BlockWeight = mxGetPr(mxCreateDoubleMatrix(winSize, winSize, mxREAL));
    double *Coor = mxGetPr(mxCreateDoubleMatrix(winSize, winSize, mxREAL)); // to label the coordinate of the pixel in the sparse matrix
    double *ColorKernel = mxGetPr(mxCreateDoubleMatrix(blockSize, blockSize, mxREAL));
    
    ARWeightJc[0] = 0;
    counter = 0;
    
    for(j=0;j<n;j++)
    {
        for(i=0;i<m;i++)
        {
            GetDepthBlock(i, j, m, n, rWin, DepthGuide, BlockDepth, Coor, InvFlag);
            
            GetColorPatch(i, j, mExtend, n, rWin, rPatch, CR, BlockR, Coor, InvFlag);
            GetColorPatch(i, j, mExtend, n, rWin, rPatch, CG, BlockG, Coor, InvFlag);
            GetColorPatch(i, j, mExtend, n, rWin, rPatch, CB, BlockB, Coor, InvFlag);

            /******* Eq.(9)******/ 
            for(q=0;q<blockSize;q++)
            {
                for(p=0;p<blockSize;p++)
                {
                    diffR = BlockR[(rWin + rPatch)*winSizeExtend + rWin + rPatch] - BlockR[(rWin + q)*winSizeExtend + rWin + p];
                    diffG = BlockG[(rWin + rPatch)*winSizeExtend + rWin + rPatch] - BlockG[(rWin + q)*winSizeExtend + rWin + p];
                    diffB = BlockB[(rWin + rPatch)*winSizeExtend + rWin + rPatch] - BlockB[(rWin + q)*winSizeExtend + rWin + p];
                    ColorKernel[q*blockSize + p] = exp(-(diffR*diffR + diffG*diffG + diffB*diffB)/(2*3*sigma4*sigma4));
                }
            }
            
            /******* Eq.(8)******/ 
            WeightSum = 0;
            
            for(t=0;t<winSize;t++)
            {
                for(s=0;s<winSize;s++)
                {
                    BlockWeight[t*winSize + s] = 0;
                    
                    if(t==rWin && s==rWin) continue;
                    
                    if(Coor[t*winSize + s]!=InvFlag)
                    {  
                        ////////// color weight /////////
                        diffColor = 0;
                        
                        for(q=0;q<blockSize;q++)
                        {
                            for(p=0;p<blockSize;p++)
                            {
                                diffR = BlockR[(rWin + q)*winSizeExtend + rWin + p] - BlockR[(t + q)*winSizeExtend + s + p];
                                diffG = BlockG[(rWin + q)*winSizeExtend + rWin + p] - BlockG[(t + q)*winSizeExtend + s + p];
                                diffB = BlockB[(rWin + q)*winSizeExtend + rWin + p] - BlockB[(t + q)*winSizeExtend + s + p];
                                
                                diffColor = diffColor + (diffR + diffG + diffB)*(diffR + diffG + diffB)*ColorKernel[q*blockSize + p];
                                
                                //diffColor = diffColor + (diffR*diffR + diffG*diffG + diffB*diffB)*ColorKernel[q*blockSize + p];
                                
                            }
                        }
                        ///// DepthWeight* ColorWeight /////
                        diffDepth = BlockDepth[t*winSize + s] - BlockDepth[rWin*winSize + rWin]; 

                        BlockWeight[t*winSize + s] = exp(-diffColor/(3*2*sigma2*sigma2))*exp(-diffDepth*diffDepth/(2*sigma1*sigma1));
                        WeightSum = BlockWeight[t*winSize + s] + WeightSum;                     
                    }
                }
            }
            
            for(t=0;t<winSize;t++)
            {
                for(s=0;s<winSize;s++)
                {
                    if(Coor[t*winSize + s]!=InvFlag)
                    {
                        if(t==rWin && s==rWin) 
                        {
                            BlockWeight[t*winSize + s] = 1;
                            
                            ARWeightPr[counter] = BlockWeight[t*winSize + s];               
                            ARWeightIr[counter] = Coor[t*winSize + s];
                            counter = counter + 1;
                            continue;
                        }
                        
                        BlockWeight[t*winSize + s] = BlockWeight[t*winSize + s]/WeightSum;
                        
                        ARWeightPr[counter] = -BlockWeight[t*winSize + s];  
                        ARWeightIr[counter] = Coor[t*winSize + s];
                        counter = counter + 1;
   
                    }
                }
            }
            
            ARWeightJc[j*m + i + 1] = counter;
        }
        
    }
}


/////////////////////////////////////////////   functions ///////////////////////////////////
void GetDepthBlock(int cur_i, int cur_j, int m, int n, int winR, double *Image, double *Block, double *Coor, int InvFlag)
{
    int i, j, sMin, sMax, tMin, tMax, winSize;
    winSize = 2*winR + 1;
    
    for(i=0;i<winSize;i++)
    {
        for(j=0;j<winSize;j++)
        {
            Block[j*winSize + i] = 0;
            Coor[j*winSize + i] = InvFlag;
        }
    }
    
    if(cur_i < winR) sMin = winR - cur_i; 
    else sMin = 0;
    
    if(cur_i + winR > m - 1) sMax = winR + m - 1 - cur_i;
    else sMax = winSize - 1;
    
    if(cur_j < winR) tMin = winR - cur_j;
    else tMin = 0;
    
    if(cur_j + winR > n - 1) tMax = winR + n - 1 - cur_j;
    else tMax = winSize - 1;
    
    for(i=sMin;i<=sMax;i++)
    {
        for(j=tMin;j<=tMax;j++)
        {
            Block[j*winSize + i] = Image[(cur_j - winR + j)*m + cur_i - winR + i];
            Coor[j*winSize + i] = (cur_j - winR + j)*m + cur_i - winR + i;
        }
    }

}

////////////////////////////////////////////////
void GetColorPatch(int cur_i, int cur_j, int mExtend, int n, int winR, int rPatch, double *Image, double *Block, double *Coor, int InvFlag)
{
    int i, j, s, t, winSize, blockSize, winSizeExtend;
    winSize = 2*winR + 1;
    blockSize = 2*rPatch + 1;
    winSizeExtend = winSize + 2*rPatch;
    
    
    for(i=0;i<winSizeExtend;i++)
    {
        for(j=0;j<winSizeExtend;j++)
        {
             Block[j*winSizeExtend + i] = 0;
        }
    }
    
    for(i=0;i<winSize;i++)
    {
        for(j=0;j<winSize;j++)
        {
            if(Coor[j*winSize + i] != InvFlag)
            {
                for(s=0;s<blockSize;s++)
                {
                    for(t=0;t<blockSize;t++)
                    {
                        Block[(j + t)*winSizeExtend + i + s] = Image[(cur_j - winR + j + t)*mExtend + cur_i - winR + i + s];
                    }
                }
            }
            
        }
    }
}