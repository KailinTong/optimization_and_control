  function targMap = targDataMap(),

  ;%***********************
  ;% Create Parameter Map *
  ;%***********************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 6;
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
    ;% Auto data (MPC_P)
    ;%
      section.nData     = 24;
      section.data(24)  = dumData; %prealloc
      
	  ;% MPC_P.A_iq
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% MPC_P.Fx
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 445;
	
	  ;% MPC_P.Fy
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 525;
	
	  ;% MPC_P.Gx
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 565;
	
	  ;% MPC_P.Gy
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 645;
	
	  ;% MPC_P.Hy
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 685;
	
	  ;% MPC_P.L
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 765;
	
	  ;% MPC_P.Qy
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 773;
	
	  ;% MPC_P.W_delta
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 1173;
	
	  ;% MPC_P.gainh1
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 1198;
	
	  ;% MPC_P.gainh2
	  section.data(11).logicalSrcIdx = 10;
	  section.data(11).dtTransOffset = 1199;
	
	  ;% MPC_P.h_02
	  section.data(12).logicalSrcIdx = 11;
	  section.data(12).dtTransOffset = 1200;
	
	  ;% MPC_P.offseth1
	  section.data(13).logicalSrcIdx = 12;
	  section.data(13).dtTransOffset = 1201;
	
	  ;% MPC_P.offseth2
	  section.data(14).logicalSrcIdx = 13;
	  section.data(14).dtTransOffset = 1202;
	
	  ;% MPC_P.roh
	  section.data(15).logicalSrcIdx = 14;
	  section.data(15).dtTransOffset = 1203;
	
	  ;% MPC_P.u_max
	  section.data(16).logicalSrcIdx = 15;
	  section.data(16).dtTransOffset = 1204;
	
	  ;% MPC_P.u_min
	  section.data(17).logicalSrcIdx = 16;
	  section.data(17).dtTransOffset = 1208;
	
	  ;% MPC_P.ue
	  section.data(18).logicalSrcIdx = 17;
	  section.data(18).dtTransOffset = 1212;
	
	  ;% MPC_P.x_max
	  section.data(19).logicalSrcIdx = 18;
	  section.data(19).dtTransOffset = 1214;
	
	  ;% MPC_P.x_min
	  section.data(20).logicalSrcIdx = 19;
	  section.data(20).dtTransOffset = 1254;
	
	  ;% MPC_P.xe
	  section.data(21).logicalSrcIdx = 20;
	  section.data(21).dtTransOffset = 1294;
	
	  ;% MPC_P.xe2
	  section.data(22).logicalSrcIdx = 21;
	  section.data(22).dtTransOffset = 1296;
	
	  ;% MPC_P.CompareToConstant2_const
	  section.data(23).logicalSrcIdx = 22;
	  section.data(23).dtTransOffset = 1297;
	
	  ;% MPC_P.CompareToConstant1_const
	  section.data(24).logicalSrcIdx = 23;
	  section.data(24).dtTransOffset = 1298;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(1) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% MPC_P.HILReadAnalog1_channels
	  section.data(1).logicalSrcIdx = 24;
	  section.data(1).dtTransOffset = 0;
	
	  ;% MPC_P.HILWriteAnalog2_channels
	  section.data(2).logicalSrcIdx = 25;
	  section.data(2).dtTransOffset = 4;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(2) = section;
      clear section
      
      section.nData     = 28;
      section.data(28)  = dumData; %prealloc
      
	  ;% MPC_P.HILInitialize_OOTerminate
	  section.data(1).logicalSrcIdx = 26;
	  section.data(1).dtTransOffset = 0;
	
	  ;% MPC_P.HILInitialize_OOExit
	  section.data(2).logicalSrcIdx = 27;
	  section.data(2).dtTransOffset = 1;
	
	  ;% MPC_P.HILInitialize_OOStart
	  section.data(3).logicalSrcIdx = 28;
	  section.data(3).dtTransOffset = 2;
	
	  ;% MPC_P.HILInitialize_OOEnter
	  section.data(4).logicalSrcIdx = 29;
	  section.data(4).dtTransOffset = 3;
	
	  ;% MPC_P.HILInitialize_AOFinal
	  section.data(5).logicalSrcIdx = 30;
	  section.data(5).dtTransOffset = 4;
	
	  ;% MPC_P.HILInitialize_POFinal
	  section.data(6).logicalSrcIdx = 31;
	  section.data(6).dtTransOffset = 5;
	
	  ;% MPC_P.HILInitialize_AIHigh
	  section.data(7).logicalSrcIdx = 32;
	  section.data(7).dtTransOffset = 6;
	
	  ;% MPC_P.HILInitialize_AILow
	  section.data(8).logicalSrcIdx = 33;
	  section.data(8).dtTransOffset = 7;
	
	  ;% MPC_P.HILInitialize_AOHigh
	  section.data(9).logicalSrcIdx = 34;
	  section.data(9).dtTransOffset = 8;
	
	  ;% MPC_P.HILInitialize_AOLow
	  section.data(10).logicalSrcIdx = 35;
	  section.data(10).dtTransOffset = 9;
	
	  ;% MPC_P.HILInitialize_AOInitial
	  section.data(11).logicalSrcIdx = 36;
	  section.data(11).dtTransOffset = 10;
	
	  ;% MPC_P.HILInitialize_AOWatchdog
	  section.data(12).logicalSrcIdx = 37;
	  section.data(12).dtTransOffset = 11;
	
	  ;% MPC_P.HILInitialize_POFrequency
	  section.data(13).logicalSrcIdx = 38;
	  section.data(13).dtTransOffset = 12;
	
	  ;% MPC_P.HILInitialize_POLeading
	  section.data(14).logicalSrcIdx = 39;
	  section.data(14).dtTransOffset = 13;
	
	  ;% MPC_P.HILInitialize_POTrailing
	  section.data(15).logicalSrcIdx = 40;
	  section.data(15).dtTransOffset = 14;
	
	  ;% MPC_P.HILInitialize_POInitial
	  section.data(16).logicalSrcIdx = 41;
	  section.data(16).dtTransOffset = 15;
	
	  ;% MPC_P.HILInitialize_POWatchdog
	  section.data(17).logicalSrcIdx = 42;
	  section.data(17).dtTransOffset = 16;
	
	  ;% MPC_P.Constant6_Value
	  section.data(18).logicalSrcIdx = 43;
	  section.data(18).dtTransOffset = 17;
	
	  ;% MPC_P.Constant7_Value
	  section.data(19).logicalSrcIdx = 44;
	  section.data(19).dtTransOffset = 37;
	
	  ;% MPC_P.Delay_InitialCondition
	  section.data(20).logicalSrcIdx = 45;
	  section.data(20).dtTransOffset = 39;
	
	  ;% MPC_P.DiscreteTransferFcn_NumCoef
	  section.data(21).logicalSrcIdx = 46;
	  section.data(21).dtTransOffset = 41;
	
	  ;% MPC_P.DiscreteTransferFcn_DenCoef
	  section.data(22).logicalSrcIdx = 47;
	  section.data(22).dtTransOffset = 43;
	
	  ;% MPC_P.DiscreteTransferFcn_InitialStat
	  section.data(23).logicalSrcIdx = 48;
	  section.data(23).dtTransOffset = 45;
	
	  ;% MPC_P.Saturation_UpperSat
	  section.data(24).logicalSrcIdx = 49;
	  section.data(24).dtTransOffset = 47;
	
	  ;% MPC_P.Saturation_LowerSat
	  section.data(25).logicalSrcIdx = 50;
	  section.data(25).dtTransOffset = 48;
	
	  ;% MPC_P.Constant_Value
	  section.data(26).logicalSrcIdx = 51;
	  section.data(26).dtTransOffset = 49;
	
	  ;% MPC_P.Saturation1_UpperSat
	  section.data(27).logicalSrcIdx = 52;
	  section.data(27).dtTransOffset = 50;
	
	  ;% MPC_P.Saturation1_LowerSat
	  section.data(28).logicalSrcIdx = 53;
	  section.data(28).dtTransOffset = 51;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(3) = section;
      clear section
      
      section.nData     = 7;
      section.data(7)  = dumData; %prealloc
      
	  ;% MPC_P.HILInitialize_CKChannels
	  section.data(1).logicalSrcIdx = 54;
	  section.data(1).dtTransOffset = 0;
	
	  ;% MPC_P.HILInitialize_DOWatchdog
	  section.data(2).logicalSrcIdx = 55;
	  section.data(2).dtTransOffset = 3;
	
	  ;% MPC_P.HILInitialize_EIInitial
	  section.data(3).logicalSrcIdx = 56;
	  section.data(3).dtTransOffset = 4;
	
	  ;% MPC_P.HILInitialize_POModes
	  section.data(4).logicalSrcIdx = 57;
	  section.data(4).dtTransOffset = 5;
	
	  ;% MPC_P.HILInitialize_POConfiguration
	  section.data(5).logicalSrcIdx = 58;
	  section.data(5).dtTransOffset = 6;
	
	  ;% MPC_P.HILInitialize_POAlignment
	  section.data(6).logicalSrcIdx = 59;
	  section.data(6).dtTransOffset = 7;
	
	  ;% MPC_P.HILInitialize_POPolarity
	  section.data(7).logicalSrcIdx = 60;
	  section.data(7).dtTransOffset = 8;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(4) = section;
      clear section
      
      section.nData     = 6;
      section.data(6)  = dumData; %prealloc
      
	  ;% MPC_P.HILInitialize_AIChannels
	  section.data(1).logicalSrcIdx = 61;
	  section.data(1).dtTransOffset = 0;
	
	  ;% MPC_P.HILInitialize_AOChannels
	  section.data(2).logicalSrcIdx = 62;
	  section.data(2).dtTransOffset = 8;
	
	  ;% MPC_P.HILInitialize_EIChannels
	  section.data(3).logicalSrcIdx = 63;
	  section.data(3).dtTransOffset = 16;
	
	  ;% MPC_P.HILInitialize_EIQuadrature
	  section.data(4).logicalSrcIdx = 64;
	  section.data(4).dtTransOffset = 24;
	
	  ;% MPC_P.HILInitialize_POChannels
	  section.data(5).logicalSrcIdx = 65;
	  section.data(5).dtTransOffset = 25;
	
	  ;% MPC_P.Delay_DelayLength
	  section.data(6).logicalSrcIdx = 66;
	  section.data(6).dtTransOffset = 33;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(5) = section;
      clear section
      
      section.nData     = 37;
      section.data(37)  = dumData; %prealloc
      
	  ;% MPC_P.HILInitialize_Active
	  section.data(1).logicalSrcIdx = 67;
	  section.data(1).dtTransOffset = 0;
	
	  ;% MPC_P.HILInitialize_AOTerminate
	  section.data(2).logicalSrcIdx = 68;
	  section.data(2).dtTransOffset = 1;
	
	  ;% MPC_P.HILInitialize_AOExit
	  section.data(3).logicalSrcIdx = 69;
	  section.data(3).dtTransOffset = 2;
	
	  ;% MPC_P.HILInitialize_DOTerminate
	  section.data(4).logicalSrcIdx = 70;
	  section.data(4).dtTransOffset = 3;
	
	  ;% MPC_P.HILInitialize_DOExit
	  section.data(5).logicalSrcIdx = 71;
	  section.data(5).dtTransOffset = 4;
	
	  ;% MPC_P.HILInitialize_POTerminate
	  section.data(6).logicalSrcIdx = 72;
	  section.data(6).dtTransOffset = 5;
	
	  ;% MPC_P.HILInitialize_POExit
	  section.data(7).logicalSrcIdx = 73;
	  section.data(7).dtTransOffset = 6;
	
	  ;% MPC_P.HILInitialize_CKPStart
	  section.data(8).logicalSrcIdx = 74;
	  section.data(8).dtTransOffset = 7;
	
	  ;% MPC_P.HILInitialize_CKPEnter
	  section.data(9).logicalSrcIdx = 75;
	  section.data(9).dtTransOffset = 8;
	
	  ;% MPC_P.HILInitialize_CKStart
	  section.data(10).logicalSrcIdx = 76;
	  section.data(10).dtTransOffset = 9;
	
	  ;% MPC_P.HILInitialize_CKEnter
	  section.data(11).logicalSrcIdx = 77;
	  section.data(11).dtTransOffset = 10;
	
	  ;% MPC_P.HILInitialize_AIPStart
	  section.data(12).logicalSrcIdx = 78;
	  section.data(12).dtTransOffset = 11;
	
	  ;% MPC_P.HILInitialize_AIPEnter
	  section.data(13).logicalSrcIdx = 79;
	  section.data(13).dtTransOffset = 12;
	
	  ;% MPC_P.HILInitialize_AOPStart
	  section.data(14).logicalSrcIdx = 80;
	  section.data(14).dtTransOffset = 13;
	
	  ;% MPC_P.HILInitialize_AOPEnter
	  section.data(15).logicalSrcIdx = 81;
	  section.data(15).dtTransOffset = 14;
	
	  ;% MPC_P.HILInitialize_AOStart
	  section.data(16).logicalSrcIdx = 82;
	  section.data(16).dtTransOffset = 15;
	
	  ;% MPC_P.HILInitialize_AOEnter
	  section.data(17).logicalSrcIdx = 83;
	  section.data(17).dtTransOffset = 16;
	
	  ;% MPC_P.HILInitialize_AOReset
	  section.data(18).logicalSrcIdx = 84;
	  section.data(18).dtTransOffset = 17;
	
	  ;% MPC_P.HILInitialize_DOPStart
	  section.data(19).logicalSrcIdx = 85;
	  section.data(19).dtTransOffset = 18;
	
	  ;% MPC_P.HILInitialize_DOPEnter
	  section.data(20).logicalSrcIdx = 86;
	  section.data(20).dtTransOffset = 19;
	
	  ;% MPC_P.HILInitialize_DOStart
	  section.data(21).logicalSrcIdx = 87;
	  section.data(21).dtTransOffset = 20;
	
	  ;% MPC_P.HILInitialize_DOEnter
	  section.data(22).logicalSrcIdx = 88;
	  section.data(22).dtTransOffset = 21;
	
	  ;% MPC_P.HILInitialize_DOReset
	  section.data(23).logicalSrcIdx = 89;
	  section.data(23).dtTransOffset = 22;
	
	  ;% MPC_P.HILInitialize_EIPStart
	  section.data(24).logicalSrcIdx = 90;
	  section.data(24).dtTransOffset = 23;
	
	  ;% MPC_P.HILInitialize_EIPEnter
	  section.data(25).logicalSrcIdx = 91;
	  section.data(25).dtTransOffset = 24;
	
	  ;% MPC_P.HILInitialize_EIStart
	  section.data(26).logicalSrcIdx = 92;
	  section.data(26).dtTransOffset = 25;
	
	  ;% MPC_P.HILInitialize_EIEnter
	  section.data(27).logicalSrcIdx = 93;
	  section.data(27).dtTransOffset = 26;
	
	  ;% MPC_P.HILInitialize_POPStart
	  section.data(28).logicalSrcIdx = 94;
	  section.data(28).dtTransOffset = 27;
	
	  ;% MPC_P.HILInitialize_POPEnter
	  section.data(29).logicalSrcIdx = 95;
	  section.data(29).dtTransOffset = 28;
	
	  ;% MPC_P.HILInitialize_POStart
	  section.data(30).logicalSrcIdx = 96;
	  section.data(30).dtTransOffset = 29;
	
	  ;% MPC_P.HILInitialize_POEnter
	  section.data(31).logicalSrcIdx = 97;
	  section.data(31).dtTransOffset = 30;
	
	  ;% MPC_P.HILInitialize_POReset
	  section.data(32).logicalSrcIdx = 98;
	  section.data(32).dtTransOffset = 31;
	
	  ;% MPC_P.HILInitialize_OOReset
	  section.data(33).logicalSrcIdx = 99;
	  section.data(33).dtTransOffset = 32;
	
	  ;% MPC_P.HILInitialize_DOFinal
	  section.data(34).logicalSrcIdx = 100;
	  section.data(34).dtTransOffset = 33;
	
	  ;% MPC_P.HILInitialize_DOInitial
	  section.data(35).logicalSrcIdx = 101;
	  section.data(35).dtTransOffset = 34;
	
	  ;% MPC_P.HILReadAnalog1_Active
	  section.data(36).logicalSrcIdx = 102;
	  section.data(36).dtTransOffset = 35;
	
	  ;% MPC_P.HILWriteAnalog2_Active
	  section.data(37).logicalSrcIdx = 103;
	  section.data(37).dtTransOffset = 36;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(6) = section;
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
    nTotSects     = 2;
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
    ;% Auto data (MPC_B)
    ;%
      section.nData     = 4;
      section.data(4)  = dumData; %prealloc
      
	  ;% MPC_B.Sum
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% MPC_B.Sum1
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% MPC_B.u1
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 2;
	
	  ;% MPC_B.Add
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 4;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(1) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% MPC_B.iter
	  section.data(1).logicalSrcIdx = 5;
	  section.data(1).dtTransOffset = 0;
	
	  ;% MPC_B.error
	  section.data(2).logicalSrcIdx = 6;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(2) = section;
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
    nTotSects     = 7;
    sectIdxOffset = 2;
    
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
    ;% Auto data (MPC_DW)
    ;%
      section.nData     = 13;
      section.data(13)  = dumData; %prealloc
      
	  ;% MPC_DW.Delay_DSTATE
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% MPC_DW.DiscreteTransferFcn_states
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 2;
	
	  ;% MPC_DW.HILInitialize_AIMinimums
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 4;
	
	  ;% MPC_DW.HILInitialize_AIMaximums
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 12;
	
	  ;% MPC_DW.HILInitialize_AOMinimums
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 20;
	
	  ;% MPC_DW.HILInitialize_AOMaximums
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 28;
	
	  ;% MPC_DW.HILInitialize_AOVoltages
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 36;
	
	  ;% MPC_DW.HILInitialize_FilterFrequency
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 44;
	
	  ;% MPC_DW.HILInitialize_POSortedFreqs
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 52;
	
	  ;% MPC_DW.HILInitialize_POValues
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 60;
	
	  ;% MPC_DW.HILReadAnalog1_Buffer
	  section.data(11).logicalSrcIdx = 10;
	  section.data(11).dtTransOffset = 68;
	
	  ;% MPC_DW.DiscreteTransferFcn_tmp
	  section.data(12).logicalSrcIdx = 11;
	  section.data(12).dtTransOffset = 72;
	
	  ;% MPC_DW.HILWriteAnalog2_Buffer
	  section.data(13).logicalSrcIdx = 12;
	  section.data(13).dtTransOffset = 74;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(1) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% MPC_DW.HILInitialize_Card
	  section.data(1).logicalSrcIdx = 13;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(2) = section;
      clear section
      
      section.nData     = 5;
      section.data(5)  = dumData; %prealloc
      
	  ;% MPC_DW.HILReadAnalog1_PWORK
	  section.data(1).logicalSrcIdx = 14;
	  section.data(1).dtTransOffset = 0;
	
	  ;% MPC_DW.Scope_PWORK.LoggedData
	  section.data(2).logicalSrcIdx = 15;
	  section.data(2).dtTransOffset = 1;
	
	  ;% MPC_DW.Scope1_PWORK.LoggedData
	  section.data(3).logicalSrcIdx = 16;
	  section.data(3).dtTransOffset = 2;
	
	  ;% MPC_DW.Scope2_PWORK.LoggedData
	  section.data(4).logicalSrcIdx = 17;
	  section.data(4).dtTransOffset = 3;
	
	  ;% MPC_DW.HILWriteAnalog2_PWORK
	  section.data(5).logicalSrcIdx = 18;
	  section.data(5).dtTransOffset = 4;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(3) = section;
      clear section
      
      section.nData     = 6;
      section.data(6)  = dumData; %prealloc
      
	  ;% MPC_DW.HILInitialize_ClockModes
	  section.data(1).logicalSrcIdx = 19;
	  section.data(1).dtTransOffset = 0;
	
	  ;% MPC_DW.HILInitialize_QuadratureModes
	  section.data(2).logicalSrcIdx = 20;
	  section.data(2).dtTransOffset = 3;
	
	  ;% MPC_DW.HILInitialize_InitialEICounts
	  section.data(3).logicalSrcIdx = 21;
	  section.data(3).dtTransOffset = 11;
	
	  ;% MPC_DW.HILInitialize_POModeValues
	  section.data(4).logicalSrcIdx = 22;
	  section.data(4).dtTransOffset = 19;
	
	  ;% MPC_DW.HILInitialize_POAlignValues
	  section.data(5).logicalSrcIdx = 23;
	  section.data(5).dtTransOffset = 27;
	
	  ;% MPC_DW.HILInitialize_POPolarityVals
	  section.data(6).logicalSrcIdx = 24;
	  section.data(6).dtTransOffset = 35;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(4) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% MPC_DW.HILInitialize_POSortedChans
	  section.data(1).logicalSrcIdx = 25;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(5) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% MPC_DW.initialized_not_empty
	  section.data(1).logicalSrcIdx = 26;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(6) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% MPC_DW.options
	  section.data(1).logicalSrcIdx = 27;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(7) = section;
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


  targMap.checksum0 = 3842294566;
  targMap.checksum1 = 3972289132;
  targMap.checksum2 = 2765226405;
  targMap.checksum3 = 1246815104;

