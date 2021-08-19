# Load data from STREAM and initialize PYDYON with it


def fromSTREAM(so, uqh, time=0, ion='D'):
    """
    Initialize a PYDYON unknown quantity vector with data from a
    STREAMOutput object.

    :param so:   STREAMOutput object to take data from.
    :param uqh:  PYDYON.UnknownQuantityHandler object to insert data with.
    :param time: Time index to initialize data from.
    :param ion:  Name of main ion species (to use for Ti).
    """
    dct = {
        'Te': so.eqsys.T_cold[time,0],
        'Ti': so.eqsys.W_i.getTemperature(ion)[time,0],
        'Ip': so.eqsys.I_p[time,0],
        'IMK2': so.eqsys.I_wall[time,0],
        # TODO insert all other ion densities
        f'ni{ion}': so.eqsys.n_i[ion].data[time,:,0]
    }

    return uqh.setvector(dct)

