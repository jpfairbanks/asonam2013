PARSE mentions BEGIN
  FILENAME "sandymentions.csv"
  FIELDS srcUser, destUser ENDFIELDS

  MAPVTX(srcVTX,srcUser,VTYPE_USER)
  MAPVTX(destVTX,destUser,VTYPE_USER)

  EDGE(ETYPE_EMAIL_FROM,srcVTX,destVTX,0)
  EDGE(ETYPE_EMAIL_FROM,destVTX,srcVTX,0)
END(emailEvent)

