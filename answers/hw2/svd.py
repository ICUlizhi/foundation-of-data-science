import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

def frobenius_norm(matrix):
    return np.sqrt(np.sum(matrix ** 2))

def svd_reconstruction(img_array, k):
    # Normalize image data to [0, 1]
    img_array = img_array / 255.0
    reconstructed = np.zeros_like(img_array)

    channel_ratios = []

    for i in range(3):  
        U, Sigma, VT = np.linalg.svd(img_array[:,:,i], full_matrices=False)  
        SigK = np.diag(Sigma[:k])  
        reconstructed[:,:,i] = U[:, :k] @ SigK @ VT[:k, :]  

        # Compute Frobenius norms
        original_norm = frobenius_norm(img_array[:,:,i])
        reconstructed_norm = frobenius_norm(reconstructed[:,:,i])
        ratio = reconstructed_norm / original_norm
        
        # Store the channel ratio
        channel_ratios.append((i, ratio))
    
    original_norm_total = frobenius_norm(img_array)
    reconstructed_norm_total = frobenius_norm(reconstructed)
    total_percentage = reconstructed_norm_total / original_norm_total

    # Rescale reconstructed image back to [0, 255]
    return Image.fromarray(np.clip(reconstructed * 255, 0, 255).astype(np.uint8)), total_percentage, channel_ratios

# Set up and read the original image
ks = [1, 2, 4, 16]
images = []
percentage = []
all_channel_ratios = []  # Store all channel ratios for each k
file = 'image1.png'

original_image = np.array(Image.open(file).convert('RGB'))

for k in ks:
    reconstructed_image, frobenius_percentage, channel_ratios = svd_reconstruction(original_image, k)
    images.append(reconstructed_image)
    percentage.append(frobenius_percentage)
    all_channel_ratios.append(channel_ratios)  # Save channel ratios for later output

# Plot the reconstructed images
num_rows = (len(ks) + 1) // 2  
fig, axes = plt.subplots(num_rows, 2, figsize=(10, 5 * num_rows))

for idx, (img, k) in enumerate(zip(images, ks)):
    row = idx // 2 
    col = idx % 2   
    axes[row, col].imshow(np.uint8(img))
    axes[row, col].set_title(f'k={k}')
    axes[row, col].axis('off')

# Hide any unused axes
for j in range(idx + 1, num_rows * 2):
    axes[j // 2, j % 2].axis('off')

plt.tight_layout()
plt.savefig('output1.png')
plt.show()

# Print total Frobenius norm percentages
for k, p in zip(ks, percentage):
    print(f"k={k}, Frobenius范数总百分比: {p:.3%}")

# Print channel ratios for each k at the end
for k, channel_ratios in zip(ks, all_channel_ratios):
    for channel_index, ratio in channel_ratios:
        print(f"k={k}, Channel {channel_index}: Frobenius范数百分比: {ratio:.3%}")
