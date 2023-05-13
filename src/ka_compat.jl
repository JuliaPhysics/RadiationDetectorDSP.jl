# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


@static if isdefined(KernelAbstractions, :Backend)
    import KernelAbstractions.Backend as _KA_Backend
    _ka_synchronize(kernel, kernel_ret) = KernelAbstractions.synchronize(kernel.backend)
    _ka_get_backend(X) = KernelAbstractions.get_backend(X)
else
    # KernelAbstractions < v0.9:
    import KernelAbstractions.Device as _KA_Backend
    _ka_synchronize(kernel, kernel_ret) = wait(kernel_ret)
    _ka_get_backend(X) = KernelAbstractions.get_device(X)
end
