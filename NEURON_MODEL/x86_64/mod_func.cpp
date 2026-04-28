#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _ghchan_reg(void);
extern void _kca2_reg(void);
extern void _kdrRL_reg(void);
extern void _L_Ca_reg(void);
extern void _mAHP_reg(void);
extern void _na3rp_reg(void);
extern void _naps_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"ghchan.mod\"");
    fprintf(stderr, " \"kca2.mod\"");
    fprintf(stderr, " \"kdrRL.mod\"");
    fprintf(stderr, " \"L_Ca.mod\"");
    fprintf(stderr, " \"mAHP.mod\"");
    fprintf(stderr, " \"na3rp.mod\"");
    fprintf(stderr, " \"naps.mod\"");
    fprintf(stderr, "\n");
  }
  _ghchan_reg();
  _kca2_reg();
  _kdrRL_reg();
  _L_Ca_reg();
  _mAHP_reg();
  _na3rp_reg();
  _naps_reg();
}

#if defined(__cplusplus)
}
#endif
