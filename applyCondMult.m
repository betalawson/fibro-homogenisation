function problem = applyCondMult(problem, factor)
% This function simply replaces the conductivity tensors in the provided
% problem with new tensors that have been scaled by the provided factor
problem.D_tensor.D_xx = problem.D_tensor.D_xx * factor;
problem.D_tensor.D_xy = problem.D_tensor.D_xy * factor;
problem.D_tensor.D_yy = problem.D_tensor.D_yy * factor;

end

