#include <Xlib.h>

void setcol_(int *r, int *g, int *b)
{
  XColor Kleur;
  int k;

  Kleur.blue = *b;
  Kleur.red = *r;
  Kleur.green = *g;
  
  k = (int) (Kleur.pixel);
//  printf ("%d\n", ColorCode(*r*256, *g*256, *b*256));
  k = ColorCode(*r*256, *g*256, *b*256);
  SetNumColor(k);
}
