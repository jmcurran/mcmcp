extern long gaussfaze, *gausssave;
extern void initnorm(long s);
extern double fastnorm ();
extern double GScale;
#define FastGauss ((--gaussfaze) ? GScale *  gausssave[gaussfaze] : fastnorm())
