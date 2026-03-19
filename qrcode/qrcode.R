install.packages("qrencoder")

library(qrencoder)

# Generate QR code raster
#add your own link
qr <- qrencode_raster("https://github.com/abigaildeegan/BI428-Transcriptomics-Analysis")

# Save as high-resolution PNG
png("qr_ta.png", width = 2000, height = 2000, res = 300, bg = "white")

# Remove all margins
par(mar = c(0, 0, 0, 0))

# Plot clean QR code
plot(qr, axes = FALSE, xlab = "", ylab = "", main = "")

# Close the device
dev.off()
