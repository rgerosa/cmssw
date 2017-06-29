import FWCore.ParameterSet.Config as cms

db_client = cms.EDAnalyzer("SiStripCommissioningOfflineDbClient",
  # general parameters
  FilePath         = cms.untracked.string('/tmp'),
  RunNumber        = cms.untracked.uint32(0),
  UseClientFile    = cms.untracked.bool(False),
  UploadHwConfig   = cms.untracked.bool(False),
  UploadAnalyses   = cms.untracked.bool(False),
  DisableDevices   = cms.untracked.bool(False),
  SaveClientFile   = cms.untracked.bool(True),
  SummaryXmlFile   = cms.untracked.FileInPath('DQM/SiStripCommissioningClients/data/summary.xml'),
  # individual parameters
  ApvTimingParameters      = cms.PSet(
    SkipFecUpdate = cms.bool(False),  # skip upload of APV PLL settings
    SkipFedUpdate = cms.bool(False),  # skip upload of FED frame finding threshold
    TargetDelay = cms.int32(-1)       # -1: latest tick (old default), otherwise target delay for all ticks' rising edge
  ),
  CalibrationParameters    = cms.PSet(),
  DaqScopeModeParameters   = cms.PSet(),
  FastFedCablingParameters = cms.PSet(),
  FedCablingParameters     = cms.PSet(),
  FedTimingParameters      = cms.PSet(),
  FineDelayParameters      = cms.PSet(
    cosmic =  cms.bool(True)
  ),
  LatencyParamameters      = cms.PSet(
    OptimizePerPartition = cms.bool(False)
  ),
  NoiseParameters          = cms.PSet(),
  OptoScanParameters       = cms.PSet(
    TargetGain = cms.double(0.863),   # target gain (0.863 ~ 690ADC for tickmark)
    SkipGainUpdate = cms.bool(False)  # wether to keep the gain the same as already on the db
  ),
  PedestalsParameters      = cms.PSet(
    DeadStripMax        = cms.double(10),    # number times the noise spread below mean noise
    NoisyStripMin       = cms.double(10),    # number times the noise spread above mean noise
    HighThreshold       = cms.double(5),    # analysis-wide high threshold for the fed zero suppression
    LowThreshold        = cms.double(2),    # analysis-wide low threshold for the fed zero suppression
    DisableBadStrips    = cms.bool(False),  # for experts! disables bad strips on the fed level 
    AddBadStrips				= cms.bool(False), #for experts! keep and add disabled bad strips. 
    KeepsStripsDisabled = cms.bool(False)   # for experts! keep strips disabled as in the db's current state
  ),
  PedsOnlyParameters       = cms.PSet(),
  ### Bad channel analysis                           
  PedsFullNoiseParameters  = cms.PSet(
        #### selections used to define a bad strip
        MaxDriftResidualCut = cms.double(20), ## the strip baseline can drift during run .. if more then N ADC count, mark the strip as bad
        MinStripNoiseCut  = cms.double(2), ## if a strip has a noise value less the N ADC, mark as low noisy i.e. bad
        MaxStripNoiseCut  = cms.double(30), ## if a strip has a noise value larger than N ADC, mark strip has high noisy i.e. bad
        MaxStripNoiseSignificanceCut = cms.double(10), ## if a strip has a noise significance larger than N, mark it as bad
        AdProbabCut   = cms.double(0.002699796063), ## this is 3 sigma quantile selection on the AndersonDarling p-value
        KsProbabCut   = cms.double(0.002699796063), ## this is 3 sigma quantile selection on the Kolmogorov Smirnov p-value
        JbProbabCut   = cms.double(0.002699796063), ## this is 3 sigma quantile selection on the jacque-Bera p-value 
        Chi2ProbabCut = cms.double(0.002699796063), ## this is 3 sigma quantile selection on the chi2 p-value (from a Gaussian fit)
        KurtosisCut   = cms.double(2), ## max value of kurtosis to identify strips with long tails
        IntegralTailCut  = cms.double(0.000000573303), ## this is the integral of the residuals placed at 5 sigma from the residual mean value
        #### Zero suppression information
        HighThreshold    = cms.double(5),  ## analysis-wide high threshold for the fed zero suppression
        LowThreshold     = cms.double(2),  ## analysis-wide low threshold for the fed zero suppression
        #### Flags on bad strips
        KeepsStripsDisabled = cms.bool(False), ### True: if a strip is bad in the db, it will be still bad; False: if a strip was bad, now the analysis will tell us if it's bad or not
        SkipEmptyStrips =  cms.bool(True), ## in the analysis, if true strips with no data are not marked, otherwise they will be bad
        DisableBadStrips = cms.bool(True), ## when the upload is performed, strips are masked
        ),
  SamplingParameters       = cms.PSet(),
  VpspScanParameters       = cms.PSet(),
)
