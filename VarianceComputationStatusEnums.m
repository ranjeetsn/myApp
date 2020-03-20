classdef VarianceComputationStatusEnums < uint32
    
   enumeration
      VARIANCE_KNOWN                                            (1)
      VARIANCE_TO_BE_COMPUTED                                   (2)
      VARIANCE_COMPUTED_SUCCESSFUL                              (3)
      VARIANCE_COMPUTATION_MAX_ITERATIONS_HIT                   (777)
      VARIANCE_COMPUTATION_FAILED_DUE_TO_NON_IDENTIFIABILITY    (888)
      VARIANCE_COMPUTATION_FAILED_DUE_TO_UNKNOWN_REASON         (9999)
   end
   
end