import pandas as pd
import numpy as np

import xtrack as xt

def attach_beam_to_sequence(sequence, beam_to_configure=1, beam_configuration=None):
    """Attach beam to sequence

    Parameters
    ----------
    mad : pymadx.madx.Madx
        Madx object
    beam_to_configure : int
        Beam number to configure
    configuration : dict
        Beam configuration

    Returns
    -------
    None

    """
    mad = sequence._madx

    # beam energy
    mad.globals.nrj = beam_configuration['beam_energy_tot']
    particle_type = 'proton'

    if 'particle_mass' in beam_configuration.keys():
        particle_mass = beam_configuration['particle_mass']
        particle_type = 'ion'
    else:
        particle_mass = mad.globals.pmass # proton mass

    if 'particle_charge' in beam_configuration.keys():
        particle_charge = beam_configuration['particle_charge']
        particle_type = 'ion'
    else:
        particle_charge = 1.

    gamma_rel = (particle_charge*beam_configuration['beam_energy_tot'])/particle_mass
    # bv and bv_aux flags
    if beam_to_configure == 1:
        ss_beam_bv, ss_bv_aux = 1, 1
    elif beam_to_configure == 2:
        ss_beam_bv, ss_bv_aux = -1, 1
    elif beam_to_configure == 4:
        ss_beam_bv, ss_bv_aux = 1, -1
    else:
        raise ValueError("beam_to_configure must be 1, 2 or 4")

    mad.globals['bv_aux'] = ss_bv_aux
    mad.input(f'''
    beam, particle={particle_type},sequence={sequence.name},
        energy={beam_configuration['beam_energy_tot']*particle_charge},
        sigt={beam_configuration.get('beam_sigt', 0.0001)},
        bv={ss_beam_bv},
        npart={beam_configuration.get('beam_npart', 1)},
        sige={beam_configuration.get('beam_sige', 1e-6)},
        ex={beam_configuration.get('beam_norm_emit_x', 1) * 1e-6 / gamma_rel},
        ey={beam_configuration.get('beam_norm_emit_y', 1) * 1e-6 / gamma_rel},
        mass={particle_mass},
        charge={particle_charge};
    ''')

def save_lines_for_closed_orbit_reference(sequence_clockwise, sequence_anticlockwise):

    lines_co_ref = {}
    if sequence_clockwise is not None:
        name_cw = sequence_clockwise.name
        lines_co_ref[name_cw + '_co_ref'] = xt.Line.from_madx_sequence(
            sequence_clockwise,
            deferred_expressions=True,
            expressions_for_element_types=('kicker', 'hkicker', 'vkicker'),
            replace_in_expr={'bv_aux': 'bvaux_' + name_cw})

    if sequence_anticlockwise is not None:
        name_acw = sequence_anticlockwise.name
        lines_co_ref[name_acw + '_co_ref'] = xt.Line.from_madx_sequence(
            sequence_anticlockwise,
            deferred_expressions=True,
            expressions_for_element_types=('kicker', 'hkicker', 'vkicker'),
            replace_in_expr={'bv_aux': 'bvaux_' + name_acw})

    return lines_co_ref



def configure_b4_from_b2(sequence_b4, sequence_b2,
        update_globals={'bv_aux': -1, 'mylhcbeam': 4}):

    mad_b2 = sequence_b2._madx
    mad_b4 = sequence_b4._madx

    var_dicts_b2 = _get_variables_dicts(mad_b2)
    var_dicts_b4 = _get_variables_dicts(mad_b4)

    b2_const=var_dicts_b2['constants']
    b4_const=var_dicts_b4['constants']
    for nn in b2_const.keys():
        if nn[0]=='_':
            print(f'The constant {nn} cannot be assigned!')
        else:
            if nn not in b4_const.keys():
                mad_b4.input(f'const {nn}={b2_const[nn]:.50e}')

    b2_indep=var_dicts_b2['independent_variables']
    b4_indep=var_dicts_b4['independent_variables']
    for nn in b2_indep.keys():
        mad_b4.input(f'{nn}={b2_indep[nn]:.50e}')

    b2_dep=var_dicts_b2['dependent_variables_expr']
    b4_dep=var_dicts_b4['dependent_variables_expr']
    for nn in b2_dep.keys():
        mad_b4.input(f'{nn}:={str(b2_dep[nn])}')

    # bv_aux and my my lhcbeam need to be defined explicitly
    for nn in update_globals.keys():
        mad_b4.input(f'{nn}={update_globals[nn]}')

    # Attach beam
    mad_b4.use(sequence_b2.name)
    beam_command = str(sequence_b2.beam)
    assert(', bv=-1.0' in beam_command)
    beam_command = beam_command.replace(', bv=-1.0', ', bv=1.0')
    mad_b4.input(beam_command)
    mad_b4.use(sequence_b4.name)

    # CHECKS
    var_dicts_b2 = _get_variables_dicts(mad_b2)
    var_dicts_b4 = _get_variables_dicts(mad_b4)

    b2_const=var_dicts_b2['constants']
    b4_const=var_dicts_b4['constants']
    for nn in b4_const.keys():
        assert b2_const[nn] == b4_const[nn]

    for nn in b2_const.keys():
        if nn not in b4_const.keys():
            print(f'Warning: b2 const {nn}={b2_const[nn]} is not in b4.')

    b2_indep=var_dicts_b2['independent_variables']
    b4_indep=var_dicts_b4['independent_variables']
    for nn in b2_indep.keys():
        if str(nn) in list(update_globals.keys()):
            continue
        assert b4_indep[nn] == b2_indep[nn]

    for nn in b4_indep.keys():
        if nn not in b2_indep.keys():
            print(f'Warning: b4 indep {nn}={b4_indep[nn]} is not in b2.')

    b2_dep=var_dicts_b2['dependent_variables_expr']
    b4_dep=var_dicts_b4['dependent_variables_expr']
    for nn in b2_dep.keys():
        if str(nn) in list(update_globals.keys()):
            continue
        assert str(b4_dep[nn]) == str(b2_dep[nn])

    for nn in b4_dep.keys():
        if nn not in b2_dep.keys():
            print(f'Warning: b4 dep {nn}={str(b4_dep[nn])} is not in b2.')

def _get_variables_dicts(mad, expressions_as_str=True):
        variables_df = _get_variables_dataframes(mad)
        outp = {
            'constants': variables_df['constants'].to_dict()['value'],
            'independent_variables':
                variables_df['independent_variables'].to_dict()['value'],
            'dependent_variables_expr':
                variables_df['dependent_variables'].to_dict()['expression'],
            'dependent_variables_val':
                variables_df['dependent_variables'].to_dict()['value'],
            }

        outp['all_variables_val'] = {kk:outp['constants'][kk] for
                kk in outp['constants'].keys()}
        outp['all_variables_val'].update(outp['independent_variables'])
        outp['all_variables_val'].update(outp['dependent_variables_val'])


        return outp

def _get_variables_dataframes(mad, expressions_as_str=True):
    '''
    Extract the dictionary of the variables and constant pandas DF of the -X global workspace.

    Returns:
        The a dictionary containing:
        - constants_df: the pandas DF with the constants
        - independent_variables: the pandas DF with the independent variables
        - dependent_variables: the pandas DF with the dependent variables
    All the three DFs have a columns 'value' with the numerical values of the costants/variables.
    The dependent_variable_df, in addition to 'value' has the following columns:
        - 'expression': the string corrensponding to the MAD-X expression
        - 'parameters': the list of parameters used in the expression
        - 'knobs': the list of the independent variables that control
            the dependent variables. Note tha the parameters can be constants and/or dependent variables,
            whereas the 'knobs' are only independent variables.
    '''
    my_dict={}
    aux=_independent_variables_df(mad)
    import numpy as np
    independent_variables_df=aux[np.logical_not(aux['constant'])].copy()
    del independent_variables_df['constant']
    constant_df=aux[aux['constant']].copy()
    del constant_df['constant']
    my_dict['constants']=constant_df
    my_dict['independent_variables']=independent_variables_df
    my_dict['dependent_variables']=_dependent_variables_df(mad)

    if expressions_as_str:
        my_dict['dependent_variables']['expression'] = (
                    my_dict['dependent_variables']['expression'].apply(
                        str))
    return my_dict

def _independent_variables_df(mad):
    '''
    Extract the pandas DF with the independent variables of the MAD-X handle.

    Returns:
        The pandas DF of the independent variables. The columns of the DF correspond to the
        - the numerical value of the independent variable (value)
        - a boolean value to know it the variable is constant or not (constant)

    See madxp/examples/variablesExamples/000_run.py
    '''

    dep_df=_dependent_variables_df(mad)
    #if len(dep_df)>0:
    #    aux=list(dep_df['knobs'].values)
    #    aux=list(itertools.chain.from_iterable(aux))
    #fundamentalSet=set(np.unique(aux))
    independent_variable_set=set(mad.globals)-set(dep_df.index)
    my_dict={}
    for i in independent_variable_set:
        my_dict[i]={}
        if mad._libmadx.get_var_type(i)==0:
            my_dict[i]['constant']=True
            # my_dict[i]['knob']=False
        else:
            my_dict[i]['constant']=False
            # if i in fundamentalSet:
            #    my_dict[i]['knob']=True
            # else:
            #    my_dict[i]['knob']=False
        my_dict[i]['value']=mad.globals[i]

    return pd.DataFrame(my_dict).transpose()[['value','constant']].sort_index()

def _dependent_variables_df(mad):
    '''
    Extract the pandas DF with the dependent variables of the MAD-X handle.

    Returns:
        The pandas DF of the dependent variables. The columns of the DF correspond to the
        - the numerical value of the dependent variable (value)
        - the string corrensponding to the MAD-X expression (expression)
        - the list of parameters used in the expression (parameters)
        - the list of the fundamental independent variables.
            These are independent variables that control numerical values of the variable (knobs).

    See madxp/examples/variablesExamples/000_run.py
    '''
    my_dict={}
    for i in list(mad.globals):
        aux=_extract_parameters(str(mad._libmadx.get_var(i)))
        if aux!=[]:
            my_dict[i]={}
            my_dict[i]['parameters']=list(np.unique(aux))

    myhash=hash(str(my_dict))
    while True:
        for i in my_dict:
            aux=[]
            for j in my_dict[i]['parameters']:
                try:
                    aux=aux+my_dict[j]['parameters']
                except:
                    aux=aux+[j]
            my_dict[i]['knobs']=list(np.unique(aux))
        if myhash==hash(str(my_dict)):
            break
        else:
            myhash=hash(str(my_dict))

    for i in my_dict:
        for j in my_dict[i]['knobs'].copy():
            if mad._libmadx.get_var_type(j)==0:
                my_dict[i]['knobs'].remove(j)
        my_dict[i]['expression']=mad._libmadx.get_var(i)
        my_dict[i]['value']=mad.globals[i]

    if len(my_dict)>0:
        return pd.DataFrame(my_dict).transpose()[['value','expression','parameters','knobs']].sort_index()
    else:
        return pd.DataFrame()

def _extract_parameters(my_string):
    '''
    Extract all the parameters of a MAD-X expression.
    Args:
        my_string: The string of the MAD-X expression to parse.

    Returns:
        The list of the parameters present in the MAD-X expression.
    '''
    if (type(my_string)=='NoneType' or my_string==None
            or my_string=='None' or my_string=='[None]' or 'table(' in my_string):
        return []
    else:
        for i in [
        '*','->','-','/','+','^','(',')','[',']',',','\'','None']:
            my_string=my_string.replace(i,' ')
        my_list=my_string.split(' ')
        my_list=list(np.unique(my_list))
        if '' in my_list:
            my_list.remove('')
        if type(my_list)=='NoneType':
            my_list=[]
        for i in my_list.copy():
            if i.isdigit() or i[0].isdigit() or i[0]=='.':
                my_list.remove(i)
        my_list=list(set(my_list)-
        set([
            'sqrt',
            'log',
            'log10',
            'exp',
            'sin',
            'cos',
            'tan',
            'asin',
            'acos',
            'atan',
            'sinh',
            'cosh',
            'tanh',
            'sinc',
            'abs',
            'erf',
            'erfc',
            'floor',
            'ceil',
            'round',
            'frac',
            'ranf',
            'gauss',
            'tgauss']))
        return my_list