<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet id="openms-qc-stylesheet" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:template match="/">
    <html><body><h1>qcML report</h1>
      <table border="1"><tr><th>run</th></tr>
        <xsl:for-each select="//runQuality"><tr><td><xsl:value-of select="@ID"/></td></tr></xsl:for-each>
      </table>
    </body></html>
  </xsl:template>
</xsl:stylesheet>
