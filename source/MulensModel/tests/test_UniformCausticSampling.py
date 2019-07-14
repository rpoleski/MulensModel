from MulensModel import UniformCausticSampling


def test_UniformCausticSampling_simple():
    """
    Tests only very basic usage of UniformCausticSampling.
    """
    sampling = UniformCausticSampling(s=1.5, q=0.9, n_points=10)
    assert sampling.s == 1.5
    assert sampling.q == 0.9
    assert sampling.n_caustics == 1

def test_UniformCausticSampling():
    """
    Test if main functions of work correctly.
    """
    s_ = [1.1, 0.5]
    q_ = [0.1, 0.001]
    n_points = 1000
    # XXX - 3 configuration? How to pass reference values then? Maybe some
    # tests do only on last one after the loop

    for (s, q) in zip(s_, q_):
        sampling = UniformCausticSampling(s=s, q=q, n_points=n_points)
        print(sampling.get_standard_parameters(0.2, 0.5, 0., 10.))
        #print(sampling.get_x_in_x_out())
        #print(sampling.orientation_check(x_caustic_in, x_caustic_out)
        print(sampling.caustic_point(0.3))
        print(sampling.which_caustic(0.4))
        # XXX - finish above


