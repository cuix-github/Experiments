#include "rad.h"
#include <math.h>


#define kMaxPolyPoints	255
#define PI	3.1415926
#define AddVector(c,a,b) (c).x=(a).x+(b).x, (c).y=(a).y+(b).y, (c).z=(a).z+(b).z
#define SubVector(c,a,b) (c).x=(a).x-(b).x, (c).y=(a).y-(b).y, (c).z=(a).z-(b).z
#define CrossVector(c,a,b)	(c).x = (a).y*(b).z - (a).z*(b).y, \
				(c).y = (a).z*(b).x - (a).x*(b).z, \
				(c).z = (a).x*(b).y - (a).y*(b).x
#define DotVector(a,b) (a).x*(b).x + (a).y*(b).y + (a).z*(b).z
#define ScaleVector(c,s) (c).x*=(s), (c).y*=(s), (c).z*=(s)
#define NormalizeVector(n,a) 	((n)=sqrt(DotVector(a,a)), \
				(n)?((a).x/=n, (a).y/=n, (a).z/=n):0)

typedef struct {
	TView	view;	/* we only need to store one face of the hemi-cube */
	double*	topFactors;	/* delta form-factors(weight for each pixel) of the top face */
	double*	sideFactors; /* delta form-factors of the side faces */
} THemicube;

static TRadParams *params;	/* input parameters */
static THemicube hemicube;	/* one hemi-cube */
static double *formfactors;	/* a form-factor array which has the same length as the number of elements */
static double totalEnergy;	/* total emitted energy; used for convergence checking */

static const TSpectra black = { {0, 0, 0} };	/* for initialization */
static int FindShootPatch(unsigned long *shootPatch);
static void SumFactors(double* formfs, int xRes, int yRes, 
	IDENTIFIER* buf, double* deltaFactors,int startY);
static void MakeTopFactors(int hres, double* deltaFactors);
static void MakeSideFactors(int hres, double* deltaFactors);
static void ComputeFormfactors(unsigned long shootPatch);
static void DistributeRad(unsigned long shootPatch);
static void DrawHCElement(TElement* ep, IDENTIFIER color);
static void DrawViewElement(TElement* ep, TColor32b *colour);
static TColor32b SpectraToRGB(TSpectra* spectra);


/* Initialize radiosity based on the input parameters p */
void InitRad(TRadParams *p, TView **displayview, TView **hemiview)
{
	int n;
	int hRes;
	unsigned long i;
	int j;
	TPatch*	pp;
	TElement* ep;
	
	params = p;
	
	/* initialize hemi-cube */
	hemicube.view.fovx = 90;
	hemicube.view.fovy = 90;
	/* make sure hemicube resolution is an even number */
	hRes = ((int)(params->hemicubeRes/2.0+0.5))*2;
	hemicube.view.xRes = hemicube.view.yRes = hRes;
	n = hRes*hRes;
	hemicube.view.buffer = (unsigned long *) calloc(n, sizeof(IDENTIFIER));
	hemicube.view.wid=0;
	hemicube.view.near = params->worldSize*0.001;
	hemicube.view.far = params->worldSize;
	
	/* take advantage of the symmetry in the delta form-factors */
	hemicube.topFactors= (double *) calloc(n/4, sizeof(double));
	hemicube.sideFactors= (double *) calloc(n/4, sizeof(double));
	MakeTopFactors(hRes/2, hemicube.topFactors);
	MakeSideFactors(hRes/2, hemicube.sideFactors);
	
	formfactors = (double *) calloc(params->nElements, sizeof(double));
	
	/* initialize radiosity */
	pp = params->patches;
	for (i=params->nPatches; i--; pp++)
		pp->unshotRad = *(pp->emission);
	ep = params->elements;
	for (i=params->nElements; i--; ep++)
		ep->rad = *(ep->patch->emission);

	/* compute total energy */
	totalEnergy = 0;
	pp = params->patches;
	for (i=params->nPatches; i--; pp++)
		for (j=0; j<kNumberOfRadSamples; j++) 
			totalEnergy += pp->emission->samples[j] * pp->area;

            /* YIORGOS' EXTRA CODE */
        *displayview = &params->displayView;
        *hemiview = &hemicube.view;
}

/* Main iterative loop */
void DoRad()
{
	unsigned long shootPatch;
	
	while (FindShootPatch(&shootPatch)) 
	{
            printf("One\n");
		ComputeFormfactors(shootPatch);
		DistributeRad(shootPatch);
		DisplayResults(&params->displayView);
	}
	
}

/* Main iterative loop, as called from the render function*/
unsigned long shootPatch_n;

int doOneIteration(void)
/* does one radiosity iteration only, returns TRUE when finished */
{	

    if (FindShootPatch(&shootPatch_n)) 
    {      
        ComputeFormfactors(shootPatch_n);
        DistributeRad(shootPatch_n);
        return 0; /*FALSE */
    } else {
        printf("Radiosity done \n");
        return 1; /* TRUE */
    }

}

/* Clean up */
void CleanUpRad()
{
	free(hemicube.topFactors);
	free(hemicube.sideFactors);
	free(hemicube.view.buffer);
	free(formfactors);

}

/* Find the next shooting patch based on the unshot energy of each patch */
/* Return 0 if convergence is reached; otherwise, return 1 */
static int FindShootPatch(unsigned long *shootPatch)
{
	int i, j;
	double energySum, error, maxEnergySum=0;
	TPatch* ep;

	ep = params->patches;
	for (i=0; i< (int) params->nPatches; i++, ep++)
	{
		energySum =0;
		for (j=0; j<kNumberOfRadSamples; j++)
			energySum += ep->unshotRad.samples[j] * ep->area;
		
		if (energySum > maxEnergySum) 
		{
			*shootPatch = i;
			maxEnergySum = energySum;
		}
	}

	error = maxEnergySum / totalEnergy;
	/* check convergence */
	if (error < params->threshold)
		return (0);		/* converged */
	else
		return (1);
	

}

/* Find out the index to the delta form-factors array */
#define Index(i)	((i)<hres? i: (hres-1- ((i)%hres)))

/* Use the largest 32bit unsigned long for background */
/* #define kBackgroundItem 0xffffffff */
/* actually we will use only 24 bits */
#define kBackgroundItem 16777215

/* Convert a hemi-cube face to form-factors */
static void SumFactors(
double* formfs, /* output */
int xRes, int yRes, /* resolution of the hemi-cube face */
IDENTIFIER* buf, /* we only need the storage of the top hemi-cube face */
double* deltaFactors, /* delta form-factors for each hemi-cube pixel */
int startY
)
{
	int i, j;
	int ii, jj;
	IDENTIFIER *ip=buf;
	int hres = xRes/2;
        register unsigned long int current_backItem;

        if (bits_for_RGB == 24)
            current_backItem = kBackgroundItem;
        else
            current_backItem = 65535; /* 2^16 - 1 */
	for (i=startY; i<yRes; i++) 
	{
		ii= Index(i)*hres;
	  	for (j=0; j<xRes; j++, ip++) 
 			if (buf[i*xRes+j] != current_backItem)  
			{
				jj = Index(j);
				formfs[buf[i*xRes+j]] += deltaFactors[ii+jj];
			}
	}
}

/* Create the delta form-factors for the top face of hemi-cube */
/* Only need to compute 1/4 of the form-factors because of the 4-way symmetry */
static void MakeTopFactors(
int hres, /* half resolution of the face */
double* deltaFactors /* output */
)
{
    int j,k;
    double xSq , ySq, xy1Sq;
	double n= hres;
	double* wp;
	double dj, dk;
	
	wp = deltaFactors;
	for (j=0; j<hres; j++)
	{
		dj = (double)j;
		ySq = (n - (dj+0.5)) / n;
       	ySq *= ySq;
       	for ( k=0 ; k<hres ; k++ )
       	{
			dk = (double)k;
       		xSq = ( n - (dk + 0.5) ) / n;
			xSq *= xSq;
			xy1Sq =  xSq + ySq + 1.0 ;
			xy1Sq *= xy1Sq;
        	*wp++ = 1.0 / (xy1Sq * PI * n * n);
       	}
    }
}

/* Create the delta form-factors for the side face of hemi-cube */
/* Only need to compute 1/4 of the form-factors because of the 4-way symmetry */
static void MakeSideFactors(
int hres, /* half resolution of the face */
double* deltaFactors /* output */
)
{
    int j,k;
    double x, xSq , y, ySq, xy1, xy1Sq;
	double n= hres;
	double* wp;
	double dj, dk;
	
	wp = deltaFactors;
	for (j=0; j<hres; j++)
	{
		dj = (double)j;
		y = (n - (dj+0.5)) / n;
       	ySq = y*y;
       	for ( k=0 ; k<hres ; k++ )
       	{
			dk = (double)k;
       		x = ( n - (dk + 0.5) ) / n;
			xSq = x*x;
			xy1 =  xSq + ySq + 1.0 ;
			xy1Sq = xy1*xy1;
        	*wp++ = y / (xy1Sq * PI * n * n);
       	}
    }
}

/* Use drand48 instead if it is supported */
#define RandomFloat ((float)(rand())/(float)RAND_MAX)

/* Compute form-factors from the shooting patch to every elements */
static void ComputeFormfactors(unsigned long shootPatch)
{
	unsigned long i;
	TVector3f	up[5]; 
	TPoint3f	lookat[5];
	TPoint3f	center;
	TVector3f	normal, tangentU, tangentV, vec;
	int face;
	double		norm;
	TPatch*		sp;
	double*		fp;
	TElement*	ep;
        double          plane_eq[4];

	/* get the center of shootPatch */
	sp = &(params->patches[shootPatch]);
	center = sp->center;
	normal = sp->normal;
        plane_eq[0] = (double)normal.x;
        plane_eq[1] = (double)normal.y;
        plane_eq[2] = (double)normal.z;
        plane_eq[3] = (double)-(normal.x*center.x + normal.y*center.y +
                                normal.z*center.z);
	
	/* rotate the hemi-cube along the normal axis of the patch randomly */
	/* this will reduce the hemi-cube aliasing artifacts */
	do {
		vec.x = RandomFloat;
		vec.y = RandomFloat;
		vec.z = RandomFloat;
		/* get a tangent vector */
		CrossVector(tangentU, normal, vec);
		NormalizeVector(norm, tangentU);
	} while (norm==0);	/* bad choice of the random vector */
	
	/* compute tangentV */
	CrossVector(tangentV, normal, tangentU);
	
	/* assign the lookats and ups for each hemicube face */
	AddVector(lookat[0], center, normal);
	up[0] = tangentU;
	AddVector(lookat[1], center, tangentU);
	up[1] = normal;
	AddVector(lookat[2], center, tangentV);
	up[2] = normal;
	SubVector(lookat[3], center, tangentU);
	up[3] = normal;
	SubVector(lookat[4], center, tangentV);
	up[4] = normal;
	
	/* position the hemicube slightly above the center of the shooting patch */
	ScaleVector(normal, params->worldSize*0.00000001);
	AddVector(hemicube.view.camera, center, normal);
	
	/* clear the formfactors */
	fp = formfactors;
	for (i=params->nElements; i--; fp++)
		*fp = 0.0;
		
	for (face=0; face < 5; face++)
	{
		hemicube.view.lookat = lookat[face];
		hemicube.view.up = up[face];

		/* draw elements */
                if (bits_for_RGB == 24) { /* a 24-bit display */
                    BeginHCDraw(&(hemicube.view), kBackgroundItem, plane_eq);
                    for (i=0; i< params->nElements; i++)
			DrawHCElement(&params->elements[i], i);	
			/* color element i with its index */
                    EndHCDraw(&(hemicube.view));
                } else if (bits_for_RGB == 8) { /* an 8-bit display */
                        /* this is a quick hack to make it work for 8-bit
                           displays maybe a better way could be found ??? */
                    
                    part_of_id = 1; /* processing first half of polygon ids */
                    BeginHCDraw(&(hemicube.view), 255, plane_eq);
                    for (i=0; i< params->nElements; i++)
			DrawHCElement(&params->elements[i], i);	
			/* color element i with its index */
                    EndHCDraw(&(hemicube.view));
                    
                    part_of_id = 2; /* second half of polygon ids */
                    BeginHCDraw(&(hemicube.view), 255, plane_eq);
                    for (i=0; i< params->nElements; i++)
			DrawHCElement(&params->elements[i], i);	
			/* color element i with its index */
                    EndHCDraw(&(hemicube.view));
                    
                } else {
                    printf("Unexpected bits per RGB colour, exiting");
                    //exit(0);
                }
		
		/* get formfactors */
		if (face==0)
			SumFactors(formfactors, hemicube.view.xRes, hemicube.view.yRes, 
				hemicube.view.buffer, hemicube.topFactors,0);
		else
			SumFactors(formfactors, hemicube.view.xRes, hemicube.view.yRes, 
				hemicube.view.buffer, hemicube.sideFactors,
				hemicube.view.yRes/2);
	}
	
	/* compute reciprocal form-factors */
	ep = params->elements;
	fp = formfactors;
	for (i=params->nElements; i--; ep++, fp++)
	{
		*fp *= sp->area / ep->area;

		/* This is a potential source of hemi-cube aliasing */
		/* To do this right, we need to subdivide the shooting patch
		and reshoot. For now we just clip it to unity */
		if ((*fp) > 1.0) 	*fp = 1.0;	
	}

}

/* Distribute radiosity form shootPatch to every element */
/* Reset the shooter's unshot radiosity to 0 */
static void DistributeRad(unsigned long shootPatch)
{
	unsigned long i;
	int j;
	TPatch* sp;
	TElement* ep;
	double* fp;
	TSpectra deltaRad;
	double w;

	sp = &(params->patches[shootPatch]);
	
	/* distribute unshotRad to every element */
	ep = params->elements;
	fp = formfactors;
	for (i=params->nElements; i--; ep++, fp++)
	{
		if ((*fp) != 0.0) 
		{
			for (j=0; j<kNumberOfRadSamples; j++)
				 deltaRad.samples[j] = 	sp->unshotRad.samples[j] * (*fp) * 
										ep->patch->reflectance->samples[j];

			/* incremental element's radiosity and patch's unshot radiosity */
			w = ep->area/ep->patch->area;
			for (j=0; j<kNumberOfRadSamples; j++) 
			{
				ep->rad.samples[j] += deltaRad.samples[j];
				ep->patch->unshotRad.samples[j] += deltaRad.samples[j] * w;
			}
		}
	}

	/* reset shooting patch's unshot radiosity */
	sp->unshotRad = black;
}

/* Convert a TSpectra (radiosity) to a TColor32b (rgb color) */
/* Assume the first three samples of the spectra are the r, g, b colors */
/* More elaborated color space transformation could be performed here */
static TColor32b
SpectraToRGB(TSpectra* spectra)
{
	TColor32b	c;
	TSpectra	r;
	double 	max=1.0;
	int k;

	for (k=kNumberOfRadSamples; k--;) {
		if (spectra->samples[k] > max)
			max = spectra->samples[k];
	}
	/* Clip the intensity*/
	r = *spectra;
	if (max>1.0) {
		for (k=kNumberOfRadSamples; k--; )
			r.samples[k] /= max;
	}
	
	/* Convert to a 32-bit color; Assume the first 3 samples in TSpectra 
	are the r, g, b colors we want. Otherwise, do color conversion here */
	c.a= 0;
	c.r= (unsigned char) (r.samples[0] * 255.0 + 0.5);
	c.g= (unsigned char) (r.samples[1] * 255.0 + 0.5);
	c.b= (unsigned char) (r.samples[2] * 255.0 + 0.5);
	
	return c;
}

static void
GetAmbient(TSpectra* ambient)
{
	int k;
	TSpectra uSum; 

	uSum=black;
	/* Do Once: compute the average reflectance */


	/* sum (unshot radiosity * area) */
	
	/* compute ambient */
	for (k=kNumberOfRadSamples; k--; )
		ambient->samples[k] = uSum.samples[k];

}

void DisplayResults(TView* view)
{
	unsigned long i;
	register TElement* ep;
	TSpectra ambient;
	GetAmbient(&ambient);
	
	BeginViewDraw(view, 0);
	ep = params->elements;
	for (i=0; i< params->nElements; i++, ep++) {
		TColor32b	c;
		TSpectra  s;
		int k;
		/* add ambient approximation */
		if (params->addAmbient) {
			for (k=kNumberOfRadSamples; k--;)
				s.samples[k] = (ep->rad.samples[k] + (ambient.samples[k]*
					ep->patch->reflectance->samples[k]))*params->intensityScale;
		} else {
			for (k=kNumberOfRadSamples; k--; )
				s.samples[k] = ep->rad.samples[k]*params->intensityScale;
		}
		/* quantize color */
		c = SpectraToRGB(&s);
		DrawViewElement(ep, &c);
	}
			
	EndViewDraw();

}

static void
DrawHCElement(TElement* ep, IDENTIFIER color)
{
	static TPoint3f pts[kMaxPolyPoints];
	int nPts = ep->nVerts;
	int j;
	for (j=0; j<nPts; j++)
		pts[j] = params->points[ep->verts[j]];
	
	DrawHCPolygon(nPts, pts, &ep->normal, color);

}

static void
DrawViewElement(TElement* ep, TColor32b *color)
{
	static TPoint3f pts[kMaxPolyPoints];
	int nPts = ep->nVerts;
	int j;
	for (j=0; j<nPts; j++)
		pts[j] = params->points[ep->verts[j]];
	
	DrawViewPolygon(nPts, pts, &ep->normal, color);

}




