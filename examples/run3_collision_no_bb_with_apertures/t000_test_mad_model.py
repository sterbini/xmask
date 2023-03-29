from cpymad.madx import Madx

mad = Madx()

energy_gev = 6800
mylhcbeam = 1

mad.input(
f'mylhcbeam={mylhcbeam};'
'''
!! S. Fartoukh. Pedestrain, sample job for using Run III optics files
!! Updated by F.F. Van der Veken


!#######################
!# Sequence and Optics #
!#######################
option,-echo,-warn;

!System,"ln -fns /afs/cern.ch/eng/acc-models/lhc/2023 run3seq";
!System,"ln -fns /afs/cern.ch/eng/lhc/optics/runIII run3opt";
!System,"ln -fns /afs/cern.ch/work/f/fvanderv/public/MADX/layout_db/files/LHC run3aper";
!System,"ln -fns /afs/cern.ch/work/f/fvanderv/public/MADX/extra_tools extra";

REAL CONST l.TAN   = 0.0;   REAL CONST l.TANAL = l.TAN;
REAL CONST l.TANAR = l.TAN; REAL CONST l.TANC  = l.TAN;
REAL CONST l.TCT   = 1.0;   REAL CONST l.TCTH  = l.TCT;REAL CONST l.TCTVA = l.TCT;
REAL CONST l.MBAS2 = 0;     REAL CONST l.MBAW  = 0;
REAL CONST l.MBCS2 = 0;     REAL CONST l.MBLS2 = 0;
REAL CONST l.MBLW  = 0;     REAL CONST l.MBWMD = 0;
REAL CONST l.MBXWH = 0;     REAL CONST l.MBXWS = 0;
REAL CONST l.MBXWT = 0;

rematchIR7=0;
createOutput=1;
make_thin=1;


call,file="run3opt/toolkit/macro.madx";
if (mylhcbeam<4){
  call,file="run3seq/lhc.seq";
  bv_aux=1;
} else {
  call,file="run3seq/lhcb4.seq";
  bv_aux=-1;
};

if (make_thin==1){    ! Thick lattice currently does not have full aperture model
    if (mylhcbeam<4){
      call, file="run3seq/aperture/aperture_as-built.b1.madx";
      call, file="run3seq/aperture/aper_tol_as-built.b1.madx";
    };
    call, file="run3seq/aperture/aperture_as-built.b2.madx";
    call, file="run3seq/aperture/aper_tol_as-built.b2.madx";
};

''')

# mad.use('lhcb1')
# twthick_b1 = mad.twiss().dframe()

mad.input('''
   !! Slice
    Option, -echo,-warn,-info;
    slicefactor=4;
    select, flag=makethin, clear;
    select, flag=makethin, class=MB,    slice=2;
    select, flag=makethin, class=MQ,    slice=2 * slicefactor;
    select, flag=makethin, class=mqxa,  slice=32* slicefactor;
    select, flag=makethin, class=mqxb,  slice=32* slicefactor;
    select, flag=makethin, pattern=mbx\. ,   slice=4;
    select, flag=makethin, pattern=mbrb\.,   slice=4;
    select, flag=makethin, pattern=mbrc\.,   slice=4;
    select, flag=makethin, pattern=mbrs\.,   slice=4;
    select, flag=makethin, pattern=mqwa\.,   slice=4;
    select, flag=makethin, pattern=mqwb\.,   slice=4;
    select, flag=makethin, pattern=mqy\.,    slice=4* slicefactor;
    select, flag=makethin, pattern=mqm\.,    slice=4* slicefactor;
    select, flag=makethin, pattern=mqmc\.,   slice=4* slicefactor;
    select, flag=makethin, pattern=mqml\.,   slice=4* slicefactor;
    select, flag=makethin, pattern=mqtlh\.,  slice=2* slicefactor;
    select, flag=makethin, pattern=mqtli\.,  slice=2* slicefactor;
    select, flag=makethin, pattern=mqt\.  ,  slice=2* slicefactor;

    if (mylhcbeam==1){
      beam;
      use,sequence=lhcb1;
      makethin, sequence=lhcb1, makedipedge=false, style=teapot, makeendmarkers=true;
     } else {
      beam;
      use,sequence=lhcb2;
      makethin, sequence=lhcb2, makedipedge=false, style=teapot, makeendmarkers=true;
    };

    call, file="run3opt/RunIII_dev/Proton_2023/opticsfile.43";

    !! Set beam
    if (mylhcbeam<4){
    Beam,particle=proton,sequence=lhcb1,energy=NRJ,NPART=numpart,sige=esigma,sigt=zsigma, ex=exn*pmass/nrj,ey=exn*pmass/nrj;          
    };
    Beam,particle=proton,sequence=lhcb2,energy=NRJ,bv=-bv_aux,NPART=numpart,sige=esigma,sigt=zsigma, ex=exn*pmass/nrj,ey=exn*pmass/nrj;

''')

mad.use('lhcb1')
twthin_b1 = mad.twiss().dframe()