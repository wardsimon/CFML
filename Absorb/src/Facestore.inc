C
C     FACESTORE.INC 
C      Holds the common block with the temporary storage of face information
C      Kept separate as not used after setup complete
       CHARACTER*80 FACE_STORE(200)
       COMMON/FACES/FACE_STORE
