  function targMap = targDataMap(),

  ;%***********************
  ;% Create Parameter Map *
  ;%***********************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 7;
    sectIdxOffset = 0;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc paramMap
    ;%
    paramMap.nSections           = nTotSects;
    paramMap.sectIdxOffset       = sectIdxOffset;
      paramMap.sections(nTotSects) = dumSection; %prealloc
    paramMap.nTotData            = -1;
    
    ;%
    ;% Auto data (Experiment_2_moving_cart_P)
    ;%
      section.nData     = 8;
      section.data(8)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_P.gain_cart
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_P.gain_pendulum
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Experiment_2_moving_cart_P.imax
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 2;
	
	  ;% Experiment_2_moving_cart_P.ChirpSignal_T
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 3;
	
	  ;% Experiment_2_moving_cart_P.ChirpSignal_f1
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 4;
	
	  ;% Experiment_2_moving_cart_P.ChirpSignal_f2
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 5;
	
	  ;% Experiment_2_moving_cart_P.Differentiator_mu
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 6;
	
	  ;% Experiment_2_moving_cart_P.Differentiator1_mu
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 7;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(1) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_P.HILReadEncoder_channels
	  section.data(1).logicalSrcIdx = 8;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_P.HILWriteAnalog_channels
	  section.data(2).logicalSrcIdx = 9;
	  section.data(2).dtTransOffset = 2;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(2) = section;
      clear section
      
      section.nData     = 33;
      section.data(33)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_P.Gain2_Gain
	  section.data(1).logicalSrcIdx = 10;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_P.WeightedSampleTime_WtEt
	  section.data(2).logicalSrcIdx = 11;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Experiment_2_moving_cart_P.WeightedSampleTime_WtEt_b
	  section.data(3).logicalSrcIdx = 12;
	  section.data(3).dtTransOffset = 2;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_OOTerminate
	  section.data(4).logicalSrcIdx = 13;
	  section.data(4).dtTransOffset = 3;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_OOExit
	  section.data(5).logicalSrcIdx = 14;
	  section.data(5).dtTransOffset = 4;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_OOStart
	  section.data(6).logicalSrcIdx = 15;
	  section.data(6).dtTransOffset = 5;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_OOEnter
	  section.data(7).logicalSrcIdx = 16;
	  section.data(7).dtTransOffset = 6;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOFinal
	  section.data(8).logicalSrcIdx = 17;
	  section.data(8).dtTransOffset = 7;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POFinal
	  section.data(9).logicalSrcIdx = 18;
	  section.data(9).dtTransOffset = 8;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AIHigh
	  section.data(10).logicalSrcIdx = 19;
	  section.data(10).dtTransOffset = 9;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AILow
	  section.data(11).logicalSrcIdx = 20;
	  section.data(11).dtTransOffset = 10;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOHigh
	  section.data(12).logicalSrcIdx = 21;
	  section.data(12).dtTransOffset = 11;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOLow
	  section.data(13).logicalSrcIdx = 22;
	  section.data(13).dtTransOffset = 12;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOInitial
	  section.data(14).logicalSrcIdx = 23;
	  section.data(14).dtTransOffset = 13;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOWatchdog
	  section.data(15).logicalSrcIdx = 24;
	  section.data(15).dtTransOffset = 14;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POFrequency
	  section.data(16).logicalSrcIdx = 25;
	  section.data(16).dtTransOffset = 15;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POLeading
	  section.data(17).logicalSrcIdx = 26;
	  section.data(17).dtTransOffset = 16;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POTrailing
	  section.data(18).logicalSrcIdx = 27;
	  section.data(18).dtTransOffset = 17;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POInitial
	  section.data(19).logicalSrcIdx = 28;
	  section.data(19).dtTransOffset = 18;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POWatchdog
	  section.data(20).logicalSrcIdx = 29;
	  section.data(20).dtTransOffset = 19;
	
	  ;% Experiment_2_moving_cart_P.Constant_Value
	  section.data(21).logicalSrcIdx = 30;
	  section.data(21).dtTransOffset = 20;
	
	  ;% Experiment_2_moving_cart_P.TransportDelay_Delay
	  section.data(22).logicalSrcIdx = 31;
	  section.data(22).dtTransOffset = 21;
	
	  ;% Experiment_2_moving_cart_P.TransportDelay_InitOutput
	  section.data(23).logicalSrcIdx = 32;
	  section.data(23).dtTransOffset = 22;
	
	  ;% Experiment_2_moving_cart_P.Step_Time
	  section.data(24).logicalSrcIdx = 33;
	  section.data(24).dtTransOffset = 23;
	
	  ;% Experiment_2_moving_cart_P.Step_Y0
	  section.data(25).logicalSrcIdx = 34;
	  section.data(25).dtTransOffset = 24;
	
	  ;% Experiment_2_moving_cart_P.Step_YFinal
	  section.data(26).logicalSrcIdx = 35;
	  section.data(26).dtTransOffset = 25;
	
	  ;% Experiment_2_moving_cart_P.Step1_Time
	  section.data(27).logicalSrcIdx = 36;
	  section.data(27).dtTransOffset = 26;
	
	  ;% Experiment_2_moving_cart_P.Step1_Y0
	  section.data(28).logicalSrcIdx = 37;
	  section.data(28).dtTransOffset = 27;
	
	  ;% Experiment_2_moving_cart_P.Step1_YFinal
	  section.data(29).logicalSrcIdx = 38;
	  section.data(29).dtTransOffset = 28;
	
	  ;% Experiment_2_moving_cart_P.Step2_Time
	  section.data(30).logicalSrcIdx = 39;
	  section.data(30).dtTransOffset = 29;
	
	  ;% Experiment_2_moving_cart_P.Step2_Y0
	  section.data(31).logicalSrcIdx = 40;
	  section.data(31).dtTransOffset = 30;
	
	  ;% Experiment_2_moving_cart_P.Step2_YFinal
	  section.data(32).logicalSrcIdx = 41;
	  section.data(32).dtTransOffset = 31;
	
	  ;% Experiment_2_moving_cart_P.Gain_Gain
	  section.data(33).logicalSrcIdx = 42;
	  section.data(33).dtTransOffset = 32;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(3) = section;
      clear section
      
      section.nData     = 7;
      section.data(7)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_P.HILInitialize_CKChannels
	  section.data(1).logicalSrcIdx = 43;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_DOWatchdog
	  section.data(2).logicalSrcIdx = 44;
	  section.data(2).dtTransOffset = 3;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_EIInitial
	  section.data(3).logicalSrcIdx = 45;
	  section.data(3).dtTransOffset = 4;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POModes
	  section.data(4).logicalSrcIdx = 46;
	  section.data(4).dtTransOffset = 5;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POConfiguration
	  section.data(5).logicalSrcIdx = 47;
	  section.data(5).dtTransOffset = 6;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POAlignment
	  section.data(6).logicalSrcIdx = 48;
	  section.data(6).dtTransOffset = 7;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POPolarity
	  section.data(7).logicalSrcIdx = 49;
	  section.data(7).dtTransOffset = 8;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(4) = section;
      clear section
      
      section.nData     = 5;
      section.data(5)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_P.HILInitialize_AIChannels
	  section.data(1).logicalSrcIdx = 50;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOChannels
	  section.data(2).logicalSrcIdx = 51;
	  section.data(2).dtTransOffset = 8;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_EIChannels
	  section.data(3).logicalSrcIdx = 52;
	  section.data(3).dtTransOffset = 16;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_EIQuadrature
	  section.data(4).logicalSrcIdx = 53;
	  section.data(4).dtTransOffset = 24;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POChannels
	  section.data(5).logicalSrcIdx = 54;
	  section.data(5).dtTransOffset = 25;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(5) = section;
      clear section
      
      section.nData     = 37;
      section.data(37)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_P.HILInitialize_Active
	  section.data(1).logicalSrcIdx = 55;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOTerminate
	  section.data(2).logicalSrcIdx = 56;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOExit
	  section.data(3).logicalSrcIdx = 57;
	  section.data(3).dtTransOffset = 2;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_DOTerminate
	  section.data(4).logicalSrcIdx = 58;
	  section.data(4).dtTransOffset = 3;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_DOExit
	  section.data(5).logicalSrcIdx = 59;
	  section.data(5).dtTransOffset = 4;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POTerminate
	  section.data(6).logicalSrcIdx = 60;
	  section.data(6).dtTransOffset = 5;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POExit
	  section.data(7).logicalSrcIdx = 61;
	  section.data(7).dtTransOffset = 6;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_CKPStart
	  section.data(8).logicalSrcIdx = 62;
	  section.data(8).dtTransOffset = 7;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_CKPEnter
	  section.data(9).logicalSrcIdx = 63;
	  section.data(9).dtTransOffset = 8;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_CKStart
	  section.data(10).logicalSrcIdx = 64;
	  section.data(10).dtTransOffset = 9;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_CKEnter
	  section.data(11).logicalSrcIdx = 65;
	  section.data(11).dtTransOffset = 10;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AIPStart
	  section.data(12).logicalSrcIdx = 66;
	  section.data(12).dtTransOffset = 11;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AIPEnter
	  section.data(13).logicalSrcIdx = 67;
	  section.data(13).dtTransOffset = 12;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOPStart
	  section.data(14).logicalSrcIdx = 68;
	  section.data(14).dtTransOffset = 13;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOPEnter
	  section.data(15).logicalSrcIdx = 69;
	  section.data(15).dtTransOffset = 14;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOStart
	  section.data(16).logicalSrcIdx = 70;
	  section.data(16).dtTransOffset = 15;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOEnter
	  section.data(17).logicalSrcIdx = 71;
	  section.data(17).dtTransOffset = 16;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_AOReset
	  section.data(18).logicalSrcIdx = 72;
	  section.data(18).dtTransOffset = 17;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_DOPStart
	  section.data(19).logicalSrcIdx = 73;
	  section.data(19).dtTransOffset = 18;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_DOPEnter
	  section.data(20).logicalSrcIdx = 74;
	  section.data(20).dtTransOffset = 19;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_DOStart
	  section.data(21).logicalSrcIdx = 75;
	  section.data(21).dtTransOffset = 20;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_DOEnter
	  section.data(22).logicalSrcIdx = 76;
	  section.data(22).dtTransOffset = 21;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_DOReset
	  section.data(23).logicalSrcIdx = 77;
	  section.data(23).dtTransOffset = 22;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_EIPStart
	  section.data(24).logicalSrcIdx = 78;
	  section.data(24).dtTransOffset = 23;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_EIPEnter
	  section.data(25).logicalSrcIdx = 79;
	  section.data(25).dtTransOffset = 24;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_EIStart
	  section.data(26).logicalSrcIdx = 80;
	  section.data(26).dtTransOffset = 25;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_EIEnter
	  section.data(27).logicalSrcIdx = 81;
	  section.data(27).dtTransOffset = 26;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POPStart
	  section.data(28).logicalSrcIdx = 82;
	  section.data(28).dtTransOffset = 27;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POPEnter
	  section.data(29).logicalSrcIdx = 83;
	  section.data(29).dtTransOffset = 28;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POStart
	  section.data(30).logicalSrcIdx = 84;
	  section.data(30).dtTransOffset = 29;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POEnter
	  section.data(31).logicalSrcIdx = 85;
	  section.data(31).dtTransOffset = 30;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_POReset
	  section.data(32).logicalSrcIdx = 86;
	  section.data(32).dtTransOffset = 31;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_OOReset
	  section.data(33).logicalSrcIdx = 87;
	  section.data(33).dtTransOffset = 32;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_DOFinal
	  section.data(34).logicalSrcIdx = 88;
	  section.data(34).dtTransOffset = 33;
	
	  ;% Experiment_2_moving_cart_P.HILInitialize_DOInitial
	  section.data(35).logicalSrcIdx = 89;
	  section.data(35).dtTransOffset = 34;
	
	  ;% Experiment_2_moving_cart_P.HILReadEncoder_Active
	  section.data(36).logicalSrcIdx = 90;
	  section.data(36).dtTransOffset = 35;
	
	  ;% Experiment_2_moving_cart_P.HILWriteAnalog_Active
	  section.data(37).logicalSrcIdx = 91;
	  section.data(37).dtTransOffset = 36;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(6) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_P.ManualSwitch_CurrentSetting
	  section.data(1).logicalSrcIdx = 92;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(7) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (parameter)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    paramMap.nTotData = nTotData;
    


  ;%**************************
  ;% Create Block Output Map *
  ;%**************************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 1;
    sectIdxOffset = 0;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc sigMap
    ;%
    sigMap.nSections           = nTotSects;
    sigMap.sectIdxOffset       = sectIdxOffset;
      sigMap.sections(nTotSects) = dumSection; %prealloc
    sigMap.nTotData            = -1;
    
    ;%
    ;% Auto data (Experiment_2_moving_cart_B)
    ;%
      section.nData     = 12;
      section.data(12)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_B.Gain1
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_B.Sum
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Experiment_2_moving_cart_B.u
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 2;
	
	  ;% Experiment_2_moving_cart_B.TmpSignalConversionAtToWorkspac
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 3;
	
	  ;% Experiment_2_moving_cart_B.Gain
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 7;
	
	  ;% Experiment_2_moving_cart_B.iA
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 8;
	
	  ;% Experiment_2_moving_cart_B.Saturation
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 9;
	
	  ;% Experiment_2_moving_cart_B.Phi
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 10;
	
	  ;% Experiment_2_moving_cart_B.x0
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 11;
	
	  ;% Experiment_2_moving_cart_B.z
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 12;
	
	  ;% Experiment_2_moving_cart_B.x0_n
	  section.data(11).logicalSrcIdx = 10;
	  section.data(11).dtTransOffset = 13;
	
	  ;% Experiment_2_moving_cart_B.z_h
	  section.data(12).logicalSrcIdx = 11;
	  section.data(12).dtTransOffset = 14;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(1) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (signal)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    sigMap.nTotData = nTotData;
    


  ;%*******************
  ;% Create DWork Map *
  ;%*******************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 12;
    sectIdxOffset = 1;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc dworkMap
    ;%
    dworkMap.nSections           = nTotSects;
    dworkMap.sectIdxOffset       = sectIdxOffset;
      dworkMap.sections(nTotSects) = dumSection; %prealloc
    dworkMap.nTotData            = -1;
    
    ;%
    ;% Auto data (Experiment_2_moving_cart_DW)
    ;%
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_DW.Phi
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(1) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_DW.Phi_b
	  section.data(1).logicalSrcIdx = 1;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(2) = section;
      clear section
      
      section.nData     = 3;
      section.data(3)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_DW.bD
	  section.data(1).logicalSrcIdx = 2;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_DW.So_inv_lastCol
	  section.data(2).logicalSrcIdx = 3;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Experiment_2_moving_cart_DW.zp
	  section.data(3).logicalSrcIdx = 4;
	  section.data(3).dtTransOffset = 2;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(3) = section;
      clear section
      
      section.nData     = 3;
      section.data(3)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_DW.bD_e
	  section.data(1).logicalSrcIdx = 5;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_DW.So_inv_lastCol_e
	  section.data(2).logicalSrcIdx = 6;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Experiment_2_moving_cart_DW.zp_p
	  section.data(3).logicalSrcIdx = 7;
	  section.data(3).dtTransOffset = 2;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(4) = section;
      clear section
      
      section.nData     = 8;
      section.data(8)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_DW.HILInitialize_AIMinimums
	  section.data(1).logicalSrcIdx = 8;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_DW.HILInitialize_AIMaximums
	  section.data(2).logicalSrcIdx = 9;
	  section.data(2).dtTransOffset = 8;
	
	  ;% Experiment_2_moving_cart_DW.HILInitialize_AOMinimums
	  section.data(3).logicalSrcIdx = 10;
	  section.data(3).dtTransOffset = 16;
	
	  ;% Experiment_2_moving_cart_DW.HILInitialize_AOMaximums
	  section.data(4).logicalSrcIdx = 11;
	  section.data(4).dtTransOffset = 24;
	
	  ;% Experiment_2_moving_cart_DW.HILInitialize_AOVoltages
	  section.data(5).logicalSrcIdx = 12;
	  section.data(5).dtTransOffset = 32;
	
	  ;% Experiment_2_moving_cart_DW.HILInitialize_FilterFrequency
	  section.data(6).logicalSrcIdx = 13;
	  section.data(6).dtTransOffset = 40;
	
	  ;% Experiment_2_moving_cart_DW.HILInitialize_POSortedFreqs
	  section.data(7).logicalSrcIdx = 14;
	  section.data(7).dtTransOffset = 48;
	
	  ;% Experiment_2_moving_cart_DW.HILInitialize_POValues
	  section.data(8).logicalSrcIdx = 15;
	  section.data(8).dtTransOffset = 56;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(5) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_DW.HILInitialize_Card
	  section.data(1).logicalSrcIdx = 16;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(6) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_DW.TransportDelay_RWORK.modelTStart
	  section.data(1).logicalSrcIdx = 17;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(7) = section;
      clear section
      
      section.nData     = 17;
      section.data(17)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_DW.HILReadEncoder_PWORK
	  section.data(1).logicalSrcIdx = 18;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_DW.Scope2_PWORK.LoggedData
	  section.data(2).logicalSrcIdx = 19;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Experiment_2_moving_cart_DW.TransportDelay_PWORK.TUbufferPtrs
	  section.data(3).logicalSrcIdx = 20;
	  section.data(3).dtTransOffset = 5;
	
	  ;% Experiment_2_moving_cart_DW.Scope3_PWORK.LoggedData
	  section.data(4).logicalSrcIdx = 21;
	  section.data(4).dtTransOffset = 7;
	
	  ;% Experiment_2_moving_cart_DW.Scope4_PWORK.LoggedData
	  section.data(5).logicalSrcIdx = 22;
	  section.data(5).dtTransOffset = 8;
	
	  ;% Experiment_2_moving_cart_DW.ToWorkspace1_PWORK.LoggedData
	  section.data(6).logicalSrcIdx = 23;
	  section.data(6).dtTransOffset = 9;
	
	  ;% Experiment_2_moving_cart_DW.ToWorkspace3_PWORK.LoggedData
	  section.data(7).logicalSrcIdx = 24;
	  section.data(7).dtTransOffset = 10;
	
	  ;% Experiment_2_moving_cart_DW.ToWorkspace4_PWORK.LoggedData
	  section.data(8).logicalSrcIdx = 25;
	  section.data(8).dtTransOffset = 11;
	
	  ;% Experiment_2_moving_cart_DW.HILWriteAnalog_PWORK
	  section.data(9).logicalSrcIdx = 26;
	  section.data(9).dtTransOffset = 12;
	
	  ;% Experiment_2_moving_cart_DW.Scope2_PWORK_i.LoggedData
	  section.data(10).logicalSrcIdx = 27;
	  section.data(10).dtTransOffset = 13;
	
	  ;% Experiment_2_moving_cart_DW.Scope4_PWORK_b.LoggedData
	  section.data(11).logicalSrcIdx = 28;
	  section.data(11).dtTransOffset = 14;
	
	  ;% Experiment_2_moving_cart_DW.Scope8_PWORK.LoggedData
	  section.data(12).logicalSrcIdx = 29;
	  section.data(12).dtTransOffset = 15;
	
	  ;% Experiment_2_moving_cart_DW.Scope1_PWORK.LoggedData
	  section.data(13).logicalSrcIdx = 30;
	  section.data(13).dtTransOffset = 16;
	
	  ;% Experiment_2_moving_cart_DW.Scope3_PWORK_c.LoggedData
	  section.data(14).logicalSrcIdx = 31;
	  section.data(14).dtTransOffset = 17;
	
	  ;% Experiment_2_moving_cart_DW.Scope5_PWORK.LoggedData
	  section.data(15).logicalSrcIdx = 32;
	  section.data(15).dtTransOffset = 18;
	
	  ;% Experiment_2_moving_cart_DW.Scope6_PWORK.LoggedData
	  section.data(16).logicalSrcIdx = 33;
	  section.data(16).dtTransOffset = 19;
	
	  ;% Experiment_2_moving_cart_DW.Scope7_PWORK.LoggedData
	  section.data(17).logicalSrcIdx = 34;
	  section.data(17).dtTransOffset = 20;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(8) = section;
      clear section
      
      section.nData     = 7;
      section.data(7)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_DW.HILInitialize_ClockModes
	  section.data(1).logicalSrcIdx = 35;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_DW.HILInitialize_QuadratureModes
	  section.data(2).logicalSrcIdx = 36;
	  section.data(2).dtTransOffset = 3;
	
	  ;% Experiment_2_moving_cart_DW.HILInitialize_InitialEICounts
	  section.data(3).logicalSrcIdx = 37;
	  section.data(3).dtTransOffset = 11;
	
	  ;% Experiment_2_moving_cart_DW.HILInitialize_POModeValues
	  section.data(4).logicalSrcIdx = 38;
	  section.data(4).dtTransOffset = 19;
	
	  ;% Experiment_2_moving_cart_DW.HILInitialize_POAlignValues
	  section.data(5).logicalSrcIdx = 39;
	  section.data(5).dtTransOffset = 27;
	
	  ;% Experiment_2_moving_cart_DW.HILInitialize_POPolarityVals
	  section.data(6).logicalSrcIdx = 40;
	  section.data(6).dtTransOffset = 35;
	
	  ;% Experiment_2_moving_cart_DW.HILReadEncoder_Buffer
	  section.data(7).logicalSrcIdx = 41;
	  section.data(7).dtTransOffset = 43;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(9) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_DW.HILInitialize_POSortedChans
	  section.data(1).logicalSrcIdx = 42;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(10) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_DW.TransportDelay_IWORK.Tail
	  section.data(1).logicalSrcIdx = 43;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(11) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% Experiment_2_moving_cart_DW.zp_not_empty
	  section.data(1).logicalSrcIdx = 44;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Experiment_2_moving_cart_DW.zp_not_empty_n
	  section.data(2).logicalSrcIdx = 45;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(12) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (dwork)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    dworkMap.nTotData = nTotData;
    


  ;%
  ;% Add individual maps to base struct.
  ;%

  targMap.paramMap  = paramMap;    
  targMap.signalMap = sigMap;
  targMap.dworkMap  = dworkMap;
  
  ;%
  ;% Add checksums to base struct.
  ;%


  targMap.checksum0 = 1163847484;
  targMap.checksum1 = 2834647576;
  targMap.checksum2 = 2657793900;
  targMap.checksum3 = 2784945;

