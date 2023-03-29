from xmask.lhc import install_errors_placeholders_hllhc

def build_sequence(mad, mylhcbeam, **kwargs):

    # Select beam
    mad.input(f'mylhcbeam = {mylhcbeam}')

    mad.input('''
      REAL CONST l.TAN   = 0.0;   REAL CONST l.TANAL = l.TAN;
      REAL CONST l.TANAR = l.TAN; REAL CONST l.TANC  = l.TAN;
      REAL CONST l.TCT   = 1.0;   REAL CONST l.TCTH  = l.TCT;REAL CONST l.TCTVA = l.TCT;
      REAL CONST l.MBAS2 = 0;     REAL CONST l.MBAW  = 0;
      REAL CONST l.MBCS2 = 0;     REAL CONST l.MBLS2 = 0;
      REAL CONST l.MBLW  = 0;     REAL CONST l.MBWMD = 0;
      REAL CONST l.MBXWH = 0;     REAL CONST l.MBXWS = 0;
      REAL CONST l.MBXWT = 0;

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
        beam; ! needed for use
        use,sequence=lhcb1; ! seems needed before makethin
        makethin, sequence=lhcb1, makedipedge=false, style=teapot, makeendmarkers=true;
        makethin, sequence=lhcb2, makedipedge=false, style=teapot, makeendmarkers=true;
      } else {
        beam; ! needed for use
        use,sequence=lhcb1; ! seems needed before makethin
        makethin, sequence=lhcb2, makedipedge=false, style=teapot, makeendmarkers=true;
      };

      ! Add more aperture markers (moved to after slicing to avoid negative drifts)
      if (mylhcbeam<4){
        call,   file="run3aper/APERTURE_YETS_2022-2023.seq";
      } else {
        call,   file="run3aper/APERTURE_YETS_2022-2023_b4_part1.seq";
        call,   file="run3aper/APERTURE_YETS_2022-2023_b4_part2.seq";
      };

''')

def apply_optics(mad, optics_file):
    mad.call(optics_file)


def call_after_last_use(mad):
    mad.input('''
    call,file="extra/align_sepdip.madx";
    exec, align_mbxw;
    exec, align_mbrc15;
    exec, align_mbx28;
    exec, align_mbrc28;
    exec, align_mbrs;
    exec, align_mbrb;
    ''')
