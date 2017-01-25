package max_entr6;

/*********************************************************************
 *******************  Maximum entropy deconvolution  **********************
 **         (of a planar scintigram with a point spread function)
 *********************************************************************
 **
 ** This imageJ plugin is based on programs written by David Simpson in
 ** FORTRAN 77 at Southampton University in 1994.
 **
 ** Some FORTRAN comments have been retained (beginning with a "c").
 ** Notes on importing the fortran program:
 **  OX usually = OPUS output (here usually imConv)
 **  CGRAD = imConj
 **  F is usually the trial data
 **  OPUS = convolve
 **  TROPUS = correlate
 **  Arrays like XI have dimensions "switched" - in fortran they have the form
 **   array(i,k), but works better in java to have array(k,i) for the situation
 **   where i represents image pixels and k is a small integer; can then just
 **   reference a whole image with XI[k]
 ********************************************************************/
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

public class Max_Entropy implements PlugIn {
  
  ImagePlus imInput, imDefault, imPSF, imSmooth, imTrial, imConv, imConj;
  ImagePlus imTemp, imXi0, imXi1, imXi2, imSGRAD, imEta0, imEta1, imEta2;
  ImageProcessor ipInput, ipDefault, ipPSF, ipSmooth, ipTrial, ipConv, ipConj;
  ImageProcessor ipTemp, ipXi0, ipXi1, ipXi2, ipSGRAD, ipEta0, ipEta1, ipEta2;
  float[] pxInput, pxDefault, pxConv, pxConj, pxTemp, pxXi0, pxXi1, pxXi2;
  float[] pxEta0, pxEta1, pxEta2, pxSGRAD, pxTrial;
  float[][] RePSF, ImPSF;

  float ChiNow(MEParams params, float AX, int M) {

    int K, L;
    float[][] A = new float[3][3];
    float[] B = new float[3];
    float Z;

    //     Finds the combination of the search direction which minimise
    //     ALPHA*S-CHI**2 for a particular value of ALPHA and output 
    //     chi**2 at this point


    float BX = 1f - AX;
    for (K = 0; K < M; K++) {
      for (L = 0; L < M; L++) {
        A[K][L] = BX * params.C2[K][L] - AX * params.S2[K][L];
      }
      B[K] = -(BX * params.C1[K] - AX * params.S1[K]);
    }

    params.beta = ChoSol(A, B, M, params.beta);
    float W = 0;
    for (K = 0; K < M; K++) {
      Z = 0;
      for (L = 0; L < M; L++) {
        Z += params.C2[K][L] * params.beta[L];
      }
      W += params.beta[K] * (params.C1[K] + 0.5 * Z);
    }
    return 1 + W;
  }

  float[] ChoSol(float[][] A, float[] B, int N, float[] X) {

    int I, I1, J, K;

    float[][] L = new float[3][3];
    float[] BL = new float[3];
    float Z;

    //     Solves optimisation problem for particular value of the Lagrange
    //     multiplier. Apparently does this by Cholesky-Decomposition
    //     Numerical recipies has a bit on this.

    L[0][0] = sqrt(A[0][0]);
    for (I = 2; I < N + 1; I++) {
      L[I - 1][0] = A[I - 1][0] / L[0][0];
      for (J = 2; J < I + 1; J++) {
        Z = 0;
        /* The -1+1 is a reminder that I have incremented J-1: */
        for (K = 1; K < J - 1 + 1; K++) {
          Z += L[I - 1][K - 1] * L[J - 1][K - 1];
        }
        Z = A[I - 1][J - 1] - Z;
        if (J == I) {
          L[I - 1][J - 1] = sqrt(Z);
        } else {
          L[I - 1][J - 1] = Z / L[J - 1][J - 1];
        }
      }
    }

    BL[0] = B[0] / L[0][0];
    for (I = 2; I < N + 1; I++) {
      Z = 0;
      for (K = 1; K < I - 1 + 1; K++) {
        Z += L[I - 1][K - 1] * BL[K - 1];
      }
      BL[I - 1] = (B[I - 1] - Z) / L[I - 1][I - 1];
    }

    X[N - 1] = BL[N - 1] / L[N - 1][N - 1];
    for (I1 = 1; I1 < N - 1 + 1; I1++) {
      I = N - I1;
      Z = 0;
      for (K = I + 1; K < N + 1; K++) {
        Z += L[K - 1][I - 1] * X[K - 1];
      }
      X[I - 1] = (BL[I - 1] - Z) / L[I - 1][I - 1];
    }
    return X;
  }

  float[] convolve(FHT fhtIm, FHT fhtPSF) {
    // Performs a convolution with the FTs as arguments
    IJ.log("Performing opus routine...");
    FHT result = null;
    result = fhtIm.multiply(fhtPSF);
    result.inverseTransform();
    result.swapQuadrants();
    result.resetMinAndMax();
    float[] res = (float[]) result.getPixels();

    /*for (int i = 0; i < res.length; i++) {
      res[i] /= res.length;
    }*/
    return res;
  }

  float[] convolve(float[] pxImp, FHT fhtPSF) {
    // Performs a convolution of the array argument with the
    // global point spread function - the FHT of which is an arg

    // Copy the array data into a temporary image
    System.arraycopy(pxImp, 0, pxTemp, 0, pxTemp.length);

    // Now perform convolution with the temporary image
    return convolve(imTemp, fhtPSF);
  }

  float[] convolve(ImagePlus imp, FHT fhtPSF) {
    // Performs a convolution of the image argument with the
    // global point spread function - the FHT of which is an arg

    // Get FT of input image
    FHT fhtIm = null;
    ImageProcessor ipFhtIm = (ImageProcessor) imp.getProperty("FHT");
    if (ipFhtIm != null) {
      fhtIm = new FHT(ipFhtIm);
    } else {
      ImageProcessor ip = imp.getProcessor();
      fhtIm = new FHT(ip);
    }
    if (!fhtIm.powerOf2Size()) {
      return null;
    }
    if (ipFhtIm == null) {
      fhtIm.transform();
    }

    // Perform the convolution with two fourier transforms
    return convolve(fhtIm, fhtPSF);
  }

  float[] correlate(FHT fhtIm, FHT fhtPSF) {
    // Performs a convolution with the FTs as arguments
    IJ.log("Performing opus routine...");
    FHT result = null;
    result = fhtIm.conjugateMultiply(fhtPSF);
    result.inverseTransform();
    result.swapQuadrants();
    result.resetMinAndMax();
    float[] res = (float[]) result.getPixels();

    /*for (int i = 0; i < res.length; i++) {
      res[i] /= res.length;
    }*/

    return res;
  }

  float[] correlate(float[] pxImp, FHT fhtPSF) {
    // Performs a convolution of the array argument with the
    // global point spread function - the FHT of which is an arg

    // Copy the array data into a temporary image
    System.arraycopy(pxImp, 0, pxTemp, 0, pxTemp.length);

    // Now perform convolution with the temporary image
    return correlate(imTemp, fhtPSF);
  }

  float[] correlate(ImagePlus imp, FHT fhtPSF) {
    // Performs a convolution of the image argument with the
    // global point spread function - the FHT of which is an arg

    // Get FT of input image
    FHT fhtIm = null;
    ImageProcessor ipFhtIm = (ImageProcessor) imp.getProperty("FHT");
    if (ipFhtIm != null) {
      fhtIm = new FHT(ipFhtIm);
    } else {
      ImageProcessor ip = imp.getProcessor();
      fhtIm = new FHT(ip);
    }
    if (!fhtIm.powerOf2Size()) {
      return null;
    }
    if (ipFhtIm == null) {
      fhtIm.transform();
    }

    // Perform the convolution with two fourier transforms
    return correlate(fhtIm, fhtPSF);
  }

  float Dist(int M, float[][] S2, float[] beta) {

    //     Calculate length of beta using entropic metric in subspace
    float Dist = 0;
    for (int K = 0; K < M; K++) {
      float Z = 0;
      for (int L = 0; L < M; L++) {
        Z -= S2[K][L] * beta[L];
      }
      Dist += beta[K] * Z;
    }

    return Dist;
  }

  float e() {
    return (float) Math.E;
  }

  MEPrefs getPreferences(MEPrefs prefs) {

    //////////////////////////////////////////////////////////////////
    // Get the windows that could be operated upon
    //////////////////////////////////////////////////////////////////
    prefs.imIDs = WindowManager.getIDList();
    if (prefs.imIDs == null) {
      IJ.noImage();
    }

    String[] titles = new String[prefs.imIDs.length];
    for (int i = 0; i < prefs.imIDs.length; i++) {
      ImagePlus imp = WindowManager.getImage(prefs.imIDs[i]);
      if (imp != null) {
        titles[i] = imp.getTitle();
      } else {
        titles[i] = "";
      }
    }
    prefs.imTitles = titles;

    int index1 = 0, index2 = 1;
    if (index2 >= titles.length) {
      index2 = 0;
    }
    //////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////
    // Make the dialogue box
    //////////////////////////////////////////////////////////////////
    GenericDialog MEOpts = new GenericDialog("Options");

    MEOpts.addChoice("Input image: ", prefs.imTitles, prefs.imTitles[index1]);
    MEOpts.addChoice(
            "Point spread image: ",
            prefs.imTitles,
            prefs.imTitles[index2]);
    MEOpts.addNumericField("C1 (Std Dev multiplier)", prefs.sigmaC1, 2);
    MEOpts.addNumericField("C2 (Std Dev addition)", prefs.sigmaC2, 2);
    MEOpts.addCheckbox("Smooth solution as default?", true);
    MEOpts.showDialog();
    //////////////////////////////////////////////////////////////////

    prefs.wasCancelled = MEOpts.wasCanceled();
    if (!prefs.wasCancelled) {

      // Get images here, since the images are globals
      index1 = MEOpts.getNextChoiceIndex();
      index2 = MEOpts.getNextChoiceIndex();

      imInput = WindowManager.getImage(prefs.imIDs[index1]);
      ipInput = imInput.getProcessor();
      imPSF = WindowManager.getImage(prefs.imIDs[index2]);
      ipPSF = imInput.getProcessor();

      prefs.imWidth = imInput.getWidth();
      prefs.imHeight = imInput.getHeight();
      prefs.imSize = prefs.imWidth * prefs.imHeight;

      // Show image for testing
      imPSF.show();
      imPSF.setTitle("PSF");

      // Get error constants
      prefs.sigmaC1 = (float) MEOpts.getNextNumber();
      prefs.sigmaC2 = (float) MEOpts.getNextNumber();
      // Avoid div by 0 errors:
      if (prefs.sigmaC1 == 0) {
        prefs.sigmaC1 += 0.00001;
      }
      //if (prefs.sigmaC2 == 0) {
      prefs.sigmaC2 += 0.0001;
      //}

      // Check the default image setting
      prefs.smoothDefault = MEOpts.getNextBoolean();
    }
    return prefs;
  }

  FHT imFFT(ImagePlus imp) {
    //Returns the fourier (actually hartley) transform of the argument image

    // Define new frequency space transform - the method is taken
    //  straight from FFTMath.java
    FHT anFHT = null;
    ImageProcessor ipFHT = (ImageProcessor) imp.getProperty("FHT");
    if (ipFHT != null) {
      anFHT = new FHT(ipFHT);
    } else {
      ImageProcessor ip = imp.getProcessor();
      anFHT = new FHT(ip);
    }
    if (!anFHT.powerOf2Size()) {
      IJ.error("FFT Math",
              "Images must be a power of 2 size (256x256, 512x512, etc.)");
      return null;
    }
    if (ipFHT == null) {
      anFHT.transform();
    }
    return anFHT;
  }

  float[][] ImReFFT(float[] fht) {
    return ImReFFT(fht, true);
  }

  float[][] ImReFFT(float[] fht, boolean boolReal) {
    // Returns an array containing the "real" or "imaginary" FFT data from the
    //  input fast hartley transform.  Based on FHTreal in FHT.java.
    int maxN = (int) Math.pow(fht.length, 0.5);
    float[][] ReOrIm = new float[maxN][maxN];

    //Iterate rows
    for (int row = 0; row < maxN; row++) {
      int base = row * maxN;
      int offs = ((maxN - row) % maxN) * maxN;
      // Iterate columns (if outside loop for speed)
      if (boolReal) {
        // Real
        for (int c = 0; c < maxN; c++) {
          ReOrIm[row][c] = (fht[base + c] + fht[offs + ((maxN - c) % maxN)]) * 0.5f;
        }
      } else {
        // Imaginary
        for (int c = 0; c < maxN; c++) {
          ReOrIm[row][c] =
                  (-fht[base + c] + fht[offs + ((maxN - c) % maxN)]) * 0.5f;
        }
      }
    }
    return ReOrIm;
  }

  float LPDist(float X0, float Y0, float M, float C) {

    float LPDIST = (Y0 - M * X0 - C) / (sqrt(1 + sqr(M)));
    LPDIST = Math.abs(LPDIST);
    return LPDIST;
  }

  void makeDefaultImage(MEPrefs prefs, MEParams params) {

    if (prefs.smoothDefault) {
      // Make smoothed image as default solution
      // Duplicate the input image and perform a smooth...
      imDefault = new Duplicator().run(imInput);
      ipDefault = imDefault.getProcessor();
      ipDefault.smooth();
    } else {
      imDefault = NewImage.createFloatImage("Default", prefs.imWidth,
              prefs.imHeight, 1, 1);
      ipDefault = imDefault.getProcessor();
      //Set default solution (Base) and scale factor (blank)
      if (params.Blank == 0) {
        ImageStatistics stat = imDefault.getStatistics();
        params.Blank = (float) stat.mean;
      } else {
        ipDefault.min(params.Blank);
        ipDefault.max(params.Blank);
      }
    }
    pxDefault = (float[]) ipDefault.getPixels();
    // Make sure the minimum value in the base image is not too small
    ipDefault.min(0.1);
    imDefault.show();
    // Get max value to set contast
    ImageStatistics stat = imDefault.getStatistics();
    ipDefault.setMinAndMax(0, stat.max);
    imDefault.updateAndRepaintWindow();
    imDefault.setTitle("Default");
  }

  void makeFloats() {
    new ImageConverter(imInput).convertToGray32();
    ipInput = imInput.getProcessor();
    pxInput = (float[]) ipInput.getPixels();

    new ImageConverter(imPSF).convertToGray32();
    ipPSF = imPSF.getProcessor();
  }

  void maxEntropy(MEPrefs prefs, MEParams params) {
    IJ.log("Find maximum entropy solution:");

    IJ.log("1:");

    // Imported from FORTRAN
    //	Initialise Parameters for calls to subroutines
    int N = prefs.imWidth * prefs.imHeight;
    boolean SumFix = false;
    int MAX = 200;

    IJ.log("2:");
    setUpErrors(prefs, params);

    //Normalise the PSF to the total counts in the image
    ImageStatistics stat = imPSF.getStatistics();
    ipPSF.multiply( 1d /
            ((double)prefs.imHeight * (double)prefs.imWidth * stat.mean));

    IJ.log("3:");
    FHT psfFHT = imFFT(imPSF);
    // Get pixel values of the FFT
    IJ.log("3.1:");
    float[] pxFHT = (float[]) psfFHT.getPixels();
    // Get real and imaginary values
    IJ.log("3.2:");
    RePSF = ImReFFT(pxFHT);
    IJ.log("3.3:");
    ImPSF = ImReFFT(pxFHT, false);
    IJ.log("3.4:");

    IJ.log("4:");
    maxEntropyAlgorithm(prefs, params, psfFHT, SumFix, MAX);

    // Make X^2 image
    float[] chisqr = new float[N];
    for (int i=0; i<N; i++){
      //    ...first make a sign map
      if (pxTrial[i] - pxInput[i] < 0){
        chisqr[i] = -1;
      }
      else {
        chisqr[i] = 1;
      }
      //  ...now create the chi^2 data
      chisqr[i] *= sqr(pxTrial[i] - pxInput[i])/sqr(params.sigma[i]);
    }
    // Display X^2 as image
    ImageProcessor ipChiSqr = new FloatProcessor(prefs.imWidth, prefs.imWidth,
            chisqr, null);
    ImagePlus imChiSqr = new ImagePlus("ChiSquared", ipChiSqr.duplicate());
    imChiSqr.show();

    // Show errors (if any)
    showErrors(params, MAX);

    IJ.log("Done");
  }

  void maxEntropyAlgorithm(MEPrefs prefs, MEParams params,
          FHT psfFHT, boolean SumFix, int MAX) {

    IJ.log("Performing maximum entropy algorithm...");

    // Set string for logging messages and errors mainly in testing
    String message;

    // Set up final target for chi-squared (number of real data points),
    // target for first iteration, iteration number and number of search
    //  directions.
    int iteration = 0;
    int m = params.M;
    //float[] SGRAD = new float[prefs.imSize];

    // Initial trial map is default solution
    imTrial = new Duplicator().run(imDefault);
    imTrial.setTitle("Trial");
    ipTrial = imTrial.getProcessor();
    pxTrial = (float[]) ipTrial.getPixels();
    imTrial.show();

    // Image to hold convolution operations
    imConv = new Duplicator().run(imDefault);
    imConv.setTitle("Conv");
    ipConv = imConv.getProcessor();
    pxConv = (float[]) ipConv.getPixels();
    imConv.show();

    // Image to hold complex multiplications
    imConj = new Duplicator().run(imDefault);
    imConj.setTitle("Conj");
    ipConj = imConj.getProcessor();
    pxConj = (float[]) ipConj.getPixels();
    imConj.show();

    // Image to hold temporary data to prevent destruction of good data
    imTemp = new Duplicator().run(imDefault);
    imTemp.setTitle("Temp");
    ipTemp = imTemp.getProcessor();
    pxTemp = (float[]) ipTemp.getPixels();
    imTemp.show();

    // Xi images for testing
    imXi0 = new Duplicator().run(imDefault);
    imXi0.setTitle("Xi0");
    ipXi0 = imXi0.getProcessor();
    pxXi0 = (float[]) ipXi0.getPixels();
    imXi0.show();

    // Xi images for testing
    imXi1 = new Duplicator().run(imDefault);
    imXi1.setTitle("Xi1");
    ipXi1 = imXi1.getProcessor();
    pxXi1 = (float[]) ipXi1.getPixels();
    imXi1.show();

    // Xi images for testing
    imXi2 = new Duplicator().run(imDefault);
    imXi2.setTitle("Xi2");
    ipXi2 = imXi2.getProcessor();
    pxXi2 = (float[]) ipXi2.getPixels();
    imXi2.show();

    // Eta images for testing
    imEta0 = new Duplicator().run(imDefault);
    imEta0.setTitle("Eta0");
    ipEta0 = imEta0.getProcessor();
    pxEta0 = (float[]) ipEta0.getPixels();
    imEta0.show();

    // Eta images for testing
    imEta1 = new Duplicator().run(imDefault);
    imEta1.setTitle("Eta1");
    ipEta1 = imEta1.getProcessor();
    pxEta1 = (float[]) ipEta1.getPixels();
    imEta1.show();

    // Eta images for testing
    imEta2 = new Duplicator().run(imDefault);
    imEta2.setTitle("Eta2");
    ipEta2 = imEta2.getProcessor();
    pxEta2 = (float[]) ipEta2.getPixels();
    imEta2.show();

    // SGRAD images for testing
    imSGRAD = new Duplicator().run(imDefault);
    imSGRAD.setTitle("SGRAD");
    ipSGRAD = imSGRAD.getProcessor();
    pxSGRAD = (float[]) ipSGRAD.getPixels();
    imSGRAD.show();

    boolean continu = true;
    while (continu == true) {
      // Generate trial solution
      IJ.log("    1:");
      // Convolve Trial image with the PSF and put into
      // convolution image (opus)
      System.arraycopy(
              convolve(imTrial, psfFHT),
              0,
              pxConv,
              0,
              pxConv.length);

      // DES Test by removing calls to convolve and copy instead
      //System.arraycopy(pxTrial, 0, pxConv, 0, pxConv.length);

      //     Find residuals for this map,  evalute chi-squared and store
      //     contributions to chi-squared
      float A = 0;
      float sigmaSqrd = 0;
      params.CHISQ = 0;
      for (int j = 0; j < prefs.imSize; j++) {
        A = pxConv[j] - pxInput[j];
        sigmaSqrd = sqr(params.sigma[j]);
        params.CHISQ += A * A / sigmaSqrd;
        pxConv[j] = 2f * A / sigmaSqrd;
        //message = "Sigma[" + j + "] = " + params.sigma[j];
        //IJ.log(message);
      }
      // Find direction of increasing chi-squared (grad chi-squared)
      // from contrbutions to chi-squared
      IJ.log("    2:");
      //tropus(imConv, psfFHT);
      // Correlate convolved image with the PSF and put into conjugate image
      System.arraycopy(
              correlate(imConv, psfFHT),
              0,
              pxConj,
              0,
              pxConj.length);

      // DES Test by removing calls to correlate and copy instead
      //System.arraycopy(pxConv, 0, pxConj, 0, pxConj.length);

      float TEST = 0;
      float SNORM = 0;
      float CNORM = 0;
      float TNORM = 0;

      for (int i = 0; i < prefs.imSize; i++) {
      
    	// Find sum of current map
        params.XSUM += pxTrial[i];
        //   Find direction of increasing entropy (grad entropy) This
        //   will be zero when f=base
        pxSGRAD[i] = (-1f * (float) Math.log(pxTrial[i] / pxDefault[i])
                / params.Blank);
        /*   Find length of grad entropy in entropy metric.
             i.e. length of rescaled direction which allows high values to
             develop in a few iterations */
        SNORM += pxSGRAD[i] * pxSGRAD[i] * pxTrial[i];
        // Find length of grad chi-squared in entropy metric
        // CGRAD in original fortran program is output of tropus
        // which is imConj in this program
        CNORM += pxConj[i] * pxConj[i] * pxTrial[i];
        // Find scalar product of two directions in entropy metric
        TNORM += pxSGRAD[i] * pxConj[i] * pxTrial[i];
      }
      
      SNORM = sqrt(SNORM);
      CNORM = sqrt(CNORM);
      float C = 1f / CNORM;

      float B = 0f;
      if (iteration != 0) {
        //       Calculate 1 - cosine of angle between grad s and
        //   grad chi-squared. Equals zero when they are parallel.
        // Will get NaN when SNORM or CNORM = 0
        if (SNORM == 0) {
          message = "SNORM returns 0 at iteration " + iteration;
          IJ.error(message);
          //break;
        }
        if (CNORM == 0) {
          message = "CNORM returns 0 at iteration " + iteration;
          IJ.error(message);
          //break;
        }

        TEST = sqrt(0.5f * (1f - TNORM / (SNORM * CNORM)));
        A = 1 / (SNORM * 2 * TEST);
        B = 1 / (CNORM * 2 * TEST);
      } else {
        A = 1;
        B = 1 / CNORM;
      }

      //     Find first two normalised search directions.
      float[][] XI = new float[3][prefs.imSize];
      for (int i = 0; i < prefs.imSize; i++) {
        //   Grad chi-squared multiplied by map
        //   (allowing high values to develop quicker than just using grad
        //    chi**2)
        XI[0][i] = pxTrial[i] * C * pxConj[i];
        //   Map times difference between grad s and grad chi**2
        XI[1][i] = pxTrial[i] * (A * pxSGRAD[i] - B * pxConj[i]);
      }

      // Display Xi as an image for testing
      System.arraycopy(XI[0], 0, pxXi0, 0, XI[0].length);
      System.arraycopy(XI[1], 0, pxXi1, 0, XI[1].length);

      //  Calculate components of these directions which keeps sum of map = 1
      //Line 236
      if (SumFix) {
        XI = Project(0, prefs.imSize, XI);
        XI = Project(1, prefs.imSize, XI);
        // Display Xi as an image for testing
        System.arraycopy(XI[0], 0, pxXi0, 0, XI[0].length);
        System.arraycopy(XI[1], 0, pxXi1, 0, XI[1].length);
        System.arraycopy(XI[2], 0, pxXi2, 0, XI[2].length);
      }
      //     Used in calculating a quadratic approximation to the problem.
      float[][] eta = new float[3][prefs.imSize];
      //ALL OPUS(N,P,XI(1,1),ETA(1,1))
      // Convolve XI[0] with the PSF and put into eta[0]
      System.arraycopy(convolve(XI[0],psfFHT), 0, eta[0], 0, XI[0].length);
      // DES Test by removing calls to convolve and copy instead
      //System.arraycopy(XI[0], 0, eta[0], 0, XI[0].length);

      // Update eta0 image
      System.arraycopy(eta[0], 0, pxEta0, 0, pxEta0.length);
      //ALL OPUS(N,P,XI(1,2),ETA(1,2))
      // Convolve XI[1] with the PSF and put into eta[1]
      System.arraycopy(convolve(XI[1],psfFHT), 0, eta[1], 0, XI[1].length);
      // DES Test by removing calls to convolve and copy instead
      //System.arraycopy(XI[1], 0, eta[1], 0, XI[1].length);


      // Update eta1 image
      System.arraycopy(eta[1], 0, pxEta1, 0, pxEta1.length);

      //     Calculate Curvatue of Chi**2 (ie grad grad chi**2)
      for (int j = 0; j < prefs.imSize; j++) {
        //OX(J)=ETA(J,2)/(SIGMA(J)*SIGMA(J))
        pxConv[j] = eta[1][j] / sqr(params.sigma[j]);
      }
      // Line 249:
      //ALL TROPUS(N,P,OX,XI(1,3))
      // Correlate XI[2] with the PSF and put back into XI[2]
      System.arraycopy(correlate(imConv, psfFHT),0, XI[2], 0, XI[2].length);

      // DES Test by removing calls to correlate and copy instead
      //System.arraycopy(pxConv, 0, XI[2], 0, pxConv.length);

      //     Calculate third search direction.
      A = 0;
      for (int i = 0; i < prefs.imSize; i++) {
        B = pxTrial[i] * XI[2][i];
        A += B * XI[2][i];
        XI[2][i] = B;
      }
      A = (1 / sqrt(A));
      for (int i = 0; i < prefs.imSize; i++) {
        XI[2][i] *= A;
      }

      //     Calculate components of this directions which keeps sum of map = 1
      if (SumFix) {
        XI = Project(2, prefs.imSize, XI);
      }

      //     Used in calculating a quadratic approximation to the problem
      System.arraycopy(convolve(XI[2],psfFHT), 0, eta[2], 0, XI[2].length);
      // DES Test by removing calls to convolve and copy instead
      //System.arraycopy(XI[2], 0, eta[2], 0, XI[2].length);

      // Update eta2 image
      System.arraycopy(eta[2], 0, pxEta2, 0, pxEta2.length);
      //     Calculate the linear components of the quadratic approximation
      for (int k = 0; k < m; k++) {
        params.S1[k] = 0;
        params.C1[k] = 0;
        for (int i = 0; i < prefs.imSize; i++) {
          params.S1[k] += XI[k][i] * pxSGRAD[i];
          params.C1[k] += XI[k][i] * pxConj[i];
        }
        params.C1[k] = params.C1[k] / params.CHISQ;
      }

      // mealg line 282
      //     Calculate the square term components of the quadratic
      //      approximation.
      for (int k = 1; k < m + 1; k++) {
        for (int L = 1; L < k + 1; L++) {
          params.S2[k - 1][L - 1] = 0;
          params.C2[k - 1][L - 1] = 0;
          for (int i = 0; i < prefs.imSize; i++) {
            params.S2[k - 1][L - 1] -= XI[k - 1][i] * XI[L - 1][i] / pxTrial[i];
          }
          for (int j = 0; j < prefs.imSize; j++) {
            params.C2[k - 1][L - 1] += eta[k - 1][j] * eta[L - 1][j]
                    / sqr(params.sigma[j]);
          }
          params.S2[k - 1][L - 1] = params.S2[k - 1][L - 1] / params.Blank;
          params.C2[k - 1][L - 1] *= 2 / params.CHISQ;
        }
      }
      //     Fill in other components of (symmetric) matrix
      params.C2[0][1] = params.C2[1][0];
      params.C2[0][2] = params.C2[2][0];
      params.C2[1][2] = params.C2[2][1];
      params.S2[0][1] = params.S2[1][0];
      params.S2[0][2] = params.S2[2][0];
      params.S2[1][2] = params.S2[2][1];

      //     Calculate entropy.
      float S = 0f;
      for (int i = 0; i < prefs.imSize; i++) {
        S -= pxTrial[i] * Math.log(pxTrial[i] / (pxDefault[i] * e()))
                / (params.Blank * e());
      }

      // Calculate diagnostic (related how linearly dependant the sum of the
      // answer is on the sum of the data I think)
      A = S * params.Blank * e() / params.XSUM;

      //     Output diagnostics
      //if (DIAFLG > 1) {
      message = "Iteration: " + iteration
              + "  TEST: " + TEST
              + "  S: " + S;
      IJ.log(message);
      message = "CHTARG: " + params.ChTArg
              + "  CHISQ: " + params.CHISQ
              + "  XSUM:" + params.XSUM
              + "  A:" + A;
      IJ.log(message);
      message = "IERR = " + params.IERR;
      IJ.log(message);
      message = "SNORM: " + SNORM
              + "  CNORM: " + CNORM
              + " TNORM: " + TNORM;
      IJ.log(message);
      //}

      //     Starting weights of search directions.
      params.beta[0] = -0.5f * params.C1[0] / params.C2[0][0];
      params.beta[1] = 0;
      params.beta[2] = 0;

      //     Find Reasonable target value for chi-squared for this
      //     iteration then find optimal values of search directions
      //     by solving quadratic approximation to optimisation with
      //     this value of chi**2 as a constraint.  On first iteration
      //      just move in grad chi-squared direction
      if (iteration != 0) {
        params = Move3(params);
        //   Check that everything's O.K.
        if (params.IERR == 3) {
          //     Every thing isn't O.K. Move says it's not possible to get
          //     chisq=chitarg
          IJ.error("Cannot achieve Chi-squared = chitarg");
        }
      }

      A = 0;
      for (int I = 0; I < prefs.imSize; I++) {
        //	Update trial map with optimal combination of search directions
        for (int K = 0; K < params.M; K++) {
          pxTrial[I] += params.beta[K] * XI[K][I];
        }
        //	Set negative values to a small positive value
        if (pxTrial[I] < 0) {
          pxTrial[I] = 0.001f * params.Blank;
        }
        //	Calculate sum map
        A += pxTrial[I];
      }
      message = "Trial zeroth pixel value: " + pxTrial[0];
      IJ.log(message);

      // Display Xi as an image for testing
      System.arraycopy(XI[0], 0, pxXi0, 0, XI[0].length);
      System.arraycopy(XI[1], 0, pxXi1, 0, XI[1].length);
      System.arraycopy(XI[2], 0, pxXi2, 0, XI[2].length);

      //     David added this test to see if default solution fitted the data
      if ((Math.abs(params.CHISQ / params.ChiZER - 1) < 0.01) &
              iteration == 0) {
        params.IERR = 2;
        IJ.error("Default solution not fitting the data");
      }

      //     Normalise trial map if flagged
      if (SumFix) {
        for (int I = 0; I < prefs.imSize; I++) {
          pxTrial[I] = pxTrial[I] / A;
        }
      }
      iteration++;

      //     Check whether algorithm has converged (ie we've reached our
      //     target value of chi-squared and chi-squared and entropy
      //     considerations take us in opposite directions.)
      if (TEST < 0.02
              & Math.abs(params.CHISQ / params.ChiZER - 1) < 0.01) {
        //        we have a solution
        params.IERR = 0;
        continu = false;
      }
      //     Check to see if we've gone past max iterations.
      if (iteration > MAX) {
        //	 We don't have a solution
        params.IERR = 1;
        continu = false;
      }
    }
  }

  MEParams Move3(MEParams params) {

    int M = params.M;

    //     Defines extremes of ALPHA parameter. At A1 chi**2 is at minimum
    //     achievable and (I think) at A2 entropy is at the maximum achievable.
    float A1 = 0, A2 = 1, CTARG, FX;

    //     Minimum achievable value of chi**2 in the sub space (rescaled)
    float CMIN = ChiNow(params, A1, M);

    //     If minimum achievable is greater than the final target then attempt
    //     to move part of the way to the minimum achievable.
    if (CMIN * params.CHISQ > params.ChiZER) {
      CTARG = 0.5f * (1 + CMIN);
    } else {
      //     define (rescaled) target as equal to final target for chi**2
      CTARG = params.ChiZER / params.CHISQ;
    }
    //     Defines extremes of "F" parameter which will equal 0 when
    //      chi=chitarg
    float F1 = CMIN - CTARG;
    float F2 = ChiNow(params, A2, M) - CTARG;

    //     Test to check that such a crossing point exists. I.E. if both ends
    //     of the scale have f values with the same sign then there's not
    //     much in point in looking for the point where they change sign!
    if (F1 * F2 > 0) {
      params.IERR = 3;
      IJ.error("Move3 error!", "Error searching for a sign-change");
      return params;
    }

    //     FIND value of ALPHA which gives chi=chitarg
    do {
      //     Find mid point in range of ALPHA
      float ANEW = 0.5f * (A1 + A2);
      //     Find which side of new point the zero crossing point lies and
      //     restrict search to that side.
      FX = ChiNow(params, ANEW, M) - CTARG;
      if (F1 * FX > 0) {
        A1 = ANEW;
        F1 = FX;
      }
      if (F2 * FX > 0) {
        A2 = ANEW;
        F2 = FX;
      }
      //     Tests to see whether this has been achieved
    } while (Math.abs(FX) >= 0.001);

    //     Check how big a step this combination of the search directions is.
    float W = Dist(M, params.S2, params.beta);

    //     If it's too big then scale it down
    if (W > 0.1 * params.XSUM / params.Blank) {
      for (int K = 0; K < M; K++) {
        params.beta[K] *= sqrt(0.1f * params.XSUM / (params.Blank * W));
      }
    }
    params.ChTArg = CTARG * params.CHISQ;
    return params;
  }

  float[][] Project(int k, int N, float[][] XI) {

    float A = 0;
    int i = 0;
    // Alter search directions to keep sum of trial map the same
    // by chopping off component which changes to value. I.E. subtract mean
    A = 0;
    for (i = 0; i < N; i++) {
      A += XI[k][i];
    }
    A = A / N;

    for (i = 0; i < N; i++) {
      XI[k][i] -= A;
    }
    return XI;
  }

  public void run(String arg) {

    // Get options
    MEPrefs prefs = new MEPrefs();
    prefs = getPreferences(prefs);
    MEParams params = setUpParams(prefs);

    if (!prefs.wasCancelled) {

      //Make images floating point
      makeFloats();
      // Sets a flat or smoothed image as default
      makeDefaultImage(prefs, params);
      // Find maximum entropy solution
      maxEntropy(prefs, params);
    }
  }

  void setUpErrors(MEPrefs prefs, MEParams params) {
    IJ.log("Setting up errors...");

    // Mask images for testing
    ImagePlus imSigma = new Duplicator().run(imDefault);
    imSigma.setTitle("Sigma");
    ImageProcessor ipSigma = imSigma.getProcessor();
    float[] pxSigma = (float[]) ipSigma.getPixels();
    imSigma.show();

    // Set up Errors. Field of view defined by a seiers of straight lines
    // Sigma is increased for points close to the edge of the
    // FOV. This is because Max_ent has little hope of allowing for
    // activity present outside the field of view.
    for (int COL = 0; COL < prefs.imHeight; COL++) {
      for (int ROW = 0; ROW < prefs.imWidth; ROW++) {

        // Initialise edge factor to small constant (prevents divide
        // by zero error in MEALG)
        float EDGEFAC = 0.000001f;

        // Test whether out of field of view.
        // Large rectangle part of FOV Defined by lines
        // y=23, y=106 , x=10 , y=120
        float X = ROW + 1;
        float Y = COL + 1;

        boolean OUTFOV = (Y < 23);
        OUTFOV = OUTFOV | (Y >= 106);
        OUTFOV = OUTFOV | (X < 10);
        OUTFOV = OUTFOV | (X >= 120);
        // Triangular "missing chunks" defined by
        // y=-0.6111x + 40.111  ;  y=0.6111x - 39.2
        // y=-0.6111x + 168.332 ;  y=0.6111x + 88.889
        OUTFOV = OUTFOV | (Y < (-0.611 * X + 40.111));
        OUTFOV = OUTFOV | (Y < (0.611 * X - 39.2));
        OUTFOV = OUTFOV | (Y >= (-0.611 * X + 168.332));
        OUTFOV = OUTFOV | (Y >= (0.611 * X + 88.889));

        if (OUTFOV) {
          EDGEFAC = 1000;
        } else {
          // Set factors to determine the rate at data quality
          // deteriorates toward the edge of the image and width of
          // band to consider.	Currently set to be equal to parameters
          // for GE. Possibly tune for GENESIS using chi-squared images
          float A = 15f;
          float B = 2f;
          int WIDTH = 4;

          // Increase sigma towards edge of field of view.
          // Vertical and horizontal lines.
          float DIST = X - 10;
          if (DIST < WIDTH) {
            EDGEFAC += A * (WIDTH - DIST) + B;
          }

          DIST = 120 - X;
          if (DIST < WIDTH) {
            EDGEFAC += A * (WIDTH - DIST) + B;
          }

          DIST = Y - 23;
          if (DIST < WIDTH) {
            EDGEFAC += A * (WIDTH - DIST) + B;
          }

          DIST = 106 - Y;
          if (DIST < WIDTH) {
            EDGEFAC += A * (WIDTH - DIST) + B;
          }

          // Diagonal lines.

          // Find distance of point from line. Number of pixels changed by this
          // bit will be small so if speed is a problem comment this bit out
          // and maybe extend width of band for other lines.
          DIST = LPDist(X, Y, -0.611f, 40.11f);
          if (DIST < WIDTH) {
            EDGEFAC += A * (WIDTH - DIST) + B;
          }

          DIST = LPDist(X, Y, 0.611f, -39.2f);
          if (DIST < WIDTH) {
            EDGEFAC += A * (WIDTH - DIST) + B;
          }

          DIST = LPDist(X, Y, -0.611f, 168.332f);
          if (DIST < WIDTH) {
            EDGEFAC += A * (WIDTH - DIST) + B;
          }

          DIST = LPDist(X, Y, 0.611f, 88.889f);
          if (DIST < WIDTH) {
            EDGEFAC += A * (WIDTH - DIST) + B;
          }
        }

        // Calculate sigma, adding edge correction and small
        // constant to prevent divide by zero error later on.
        params.sigma[COL + ROW * prefs.imWidth] = prefs.sigmaC1
                * sqrt(pxInput[COL + ROW * prefs.imWidth])
                + prefs.sigmaC2 + EDGEFAC;
      }
    }
    System.arraycopy(params.sigma, 0, pxSigma, 0, pxSigma.length);
  }

  MEParams setUpParams(MEPrefs prefs) {
    IJ.log("Setting up parameters...");
    MEParams params = new MEParams();
    params.beta = new float[params.M]; //M defaults to 3
    params.S1 = new float[params.M];
    params.C1 = new float[params.M];
    params.sigma = new float[prefs.imSize];
    params.S2 = new float[params.M][params.M];
    params.C2 = new float[params.M][params.M];
    params.ChiZER = params.NoP;
    params.ChTArg = params.ChiZER;

    return params;
  }

  void showErrors(MEParams params, int max) {
    String message = "";
    switch (params.IERR) {
      case 1:
        message += "Did not converge in ";
        message += Integer.toString(max);
        message += " iterations. \nIncrease error constants or reduce \n";
        message += "number of smooths of data for default";
        break;
      case 2:
        message += "Default solution fitted the data. \nDecrease error";
        message += " constants or increase \nnumber of smooths of data";
        message += " for default";
        break;
      case 3:
        message += "Trial map fitted data on iteration 1.\n";
        message += "Decrease error constants or increase\n";
        message += "number of smooths of data for default";
        break;
      case 4:
        message += "Wrong number of map points N in " +
                "Max Entropy Algorithm";
        break;
      case 5:
        message += "Wrong size data array P in " +
                "Max Entropy Algorithm";
        break;
      case 6:
        message += "Wrong value of dispflg in " +
                "Max Entropy Algorithm";
        break;
      case 7:
        message += "Wrong number of data points NOP in " +
                "Max Entropy Algorithm";
        break;
    }
    if (params.IERR != 0) {
      IJ.error(message);
    }
  }

  double sqr(double number) {
    return number * number;
  }

  double[] sqr(double[] numbers) {
    for (int i = 0; i < numbers.length; i++) {
      numbers[i] = sqr(numbers[i]);
    }
    return numbers;
  }

  float sqr(float number) {
    return number * number;
  }

  float[] sqr(float[] numbers) {
    for (int i = 0; i < numbers.length; i++) {
      numbers[i] = sqr(numbers[i]);
    }
    return numbers;
  }

  float sqrt(float number) {
    return (float) Math.sqrt(number);
  }
}

/**
 * Object to hold values passed between subroutines
 * @author neil
 *
 */
class MEParams {

  public int IERR = 0, M = 3, NoP = 8650;
  public float[] beta, S1, C1, sigma;
  public float[][] S2, C2;
  public float CHISQ = 0, ChTArg, ChiZER, XSUM = 0, Blank = 1;

  public MEParams() {
  }
}

/**
 * Object to hold preferences
 */
class MEPrefs {

  public boolean wasCancelled = false;
  public boolean smoothDefault = true;
  public int[] imIDs;
  public int imHeight, imWidth, imSize;
  public String[] imTitles;
  public float sigmaC1 = 1f, sigmaC2 = 0.5f;

  public MEPrefs() {
  }
}