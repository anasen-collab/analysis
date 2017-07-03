

static Int_t mylog2 (Int_t val) {
    Int_t ret = -1;
    while (val != 0) {
        val >>= 1;
        ret++;
    }
    return ret;
}

static Int_t mylog2_low (Int_t val) {
    Int_t ret = 0;
    if (val != 0) { 
      while ((val&0x001) != 1) {
          val >>= 1;
          ret++;
      } 
    } else { 
    ret = -1; 
    }
    return ret;
}

